#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <omp.h>
#include <memory.h>
#include <math.h>

#include "graph_struct.h"

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})


void printArray(int* array, int n){
	for (int i=0; i<n; i++){
		printf("%d: %d\n", i, array[i]);
	}
}



void parallelBellmanFord(Graph* graph, int sourceVertex) {
    int V = graph->numVertices;
	int E = graph->numEdges;
	Edge** edges = graph->edges;
	
	Edge* edge;
	int src;
	int dest;
	int weight;

	int dist_src;
	int dist_dest;

	int* final_dist;
	int* final_predecessors;
	final_dist = (int*) malloc(V * sizeof(int));
	final_predecessors = (int*) malloc(V * sizeof(int));

	int num_threads = omp_get_max_threads();
	
	// creating an array of distances and predecessors per thread
	// to avoid to write on the same structure
	int **dist_per_thread;
	int **predecessors_per_thread;
	dist_per_thread = (int**) malloc(num_threads * sizeof(int*));
	predecessors_per_thread = (int**) malloc(num_threads * sizeof(int*));

	int dynamic_V = 4 * sqrt(V/num_threads);
	int dynamic_E = 4 * sqrt(E/num_threads);
	dynamic_V = max(1, dynamic_V);
	dynamic_E = max(1, dynamic_E);
	

	// Step 1: initialize graph
	#pragma omp parallel for schedule(dynamic, dynamic_V)
	for (int v = 0; v < V; v++) {
		// Initialize the distance to all vertices to infinity
		final_dist[v] = INT_MAX;
		// And having a null predecessor
		final_predecessors[v] = -1;
	}
	// The distance from the source to itself is, of course, zero
	final_dist[sourceVertex] = 0;

	// initialize the per thread arrays
	for(int t = 0; t < num_threads; t++) {
		dist_per_thread[t] = (int*) malloc(V * sizeof(int));
		predecessors_per_thread[t] = (int*) malloc(V * sizeof(int));
	}


	// Step 2: relax edges repeatedly
	for (int v = 0; v < V - 1; v++) {

		#pragma omp parallel default(none) \
		shared(V, E, edges, num_threads, dynamic_V, dynamic_E) \
		shared(final_dist, final_predecessors, dist_per_thread, predecessors_per_thread) \
		private(src, dest, weight, dist_src, dist_dest) 
		{
			int id_thread = omp_get_thread_num();

			// update the arrays of each thread
			#pragma omp for collapse(2) schedule(dynamic, dynamic_V)
			for(int t = 0; t < num_threads; t++){
				for(int v1 = 0; v1 < V; v1++) {
					dist_per_thread[t][v1] = final_dist[v1];
					predecessors_per_thread[t][v1] = final_predecessors[v1];
				}
			}
			#pragma omp barrier

			#pragma omp for schedule(dynamic, dynamic_E)
			for (int e = 0; e < E; e++) {
				src = edges[e]->src;
				dist_src = dist_per_thread[id_thread][src];

				if (dist_src != INT_MAX){
					dest = edges[e]->dest;
					weight = edges[e]->weight;
					dist_dest = dist_per_thread[id_thread][dest];

					if (dist_src + weight < dist_dest) {
						dist_per_thread[id_thread][dest] = dist_src + weight;
						predecessors_per_thread[id_thread][dest] = src;
					}
				}
			}
		
			#pragma omp barrier
			
			// aggregate the results
			#pragma omp for schedule(dynamic, dynamic_V)
			for (int v1=0; v1 < V; v1++) {
				
				for (int t=0; t < num_threads; t++) {
					if (dist_per_thread[t][v1] < final_dist[v1]){
						final_dist[v1] = dist_per_thread[t][v1];
						final_predecessors[v1] = predecessors_per_thread[t][v1];
					}
						
				}
			}
		}
	}
	
	
	// Step 3: check for negative-weight cycles
	int negative_cycle = 0;
	#pragma omp parallel for reduction(max:negative_cycle) private(src) private(dest) private(weight) schedule(dynamic, dynamic_E)
    for (int j = 0; j < E; j++) {
		Edge* edge = graph->edges[j];
		src = edge->src;
		dest = edge->dest;
		weight = edge->weight;
        
		if (final_dist[src] != INT_MAX && final_dist[src] + weight < final_dist[dest]) {
			negative_cycle = 1;  
		}
		else{
			negative_cycle = 0;
		}
    }

	if (negative_cycle > 0){
		printf("GRAPH CONTAINS A NEGATIVE-WEIGHT CYCLE\n");
		return;
	}

	free(final_dist);
	free(final_predecessors);
	free(dist_per_thread);
	free(predecessors_per_thread);

    return;  // No negative-weight cycles
}



void parallelBellmanFordMemcpy(Graph* graph, int sourceVertex, int dyn) {
    int V = graph->numVertices;
	int E = graph->numEdges;
	Edge** edges = graph->edges;
	
	Edge* edge;
	int src;
	int dest;
	int weight;

	int dist_src;
	int dist_dest;

	int* final_dist;
	int* final_predecessors;
	final_dist = (int*) malloc(V * sizeof(int));
	final_predecessors = (int*) malloc(V * sizeof(int));

	int num_threads = omp_get_max_threads();
	
	// creating an array of distances and predecessors per thread
	// to avoid to write on the same structure
	int **dist_per_thread;
	int **predecessors_per_thread;
	dist_per_thread = (int**) malloc(num_threads * sizeof(int*));
	predecessors_per_thread = (int**) malloc(num_threads * sizeof(int*));

	int dynamic_V, dynamic_E;
	if (dyn == 1){
		dynamic_V = 1;
		dynamic_E = 1;
	}
	else{
		dynamic_V = 4 * sqrt(V/num_threads);
		dynamic_E = 4 * sqrt(E/num_threads);
		dynamic_V = max(1, dynamic_V);
		dynamic_E = max(1, dynamic_E);
	}

	
	// Step 1: initialize graph
	#pragma omp parallel for schedule(dynamic, dynamic_V)
	for (int v = 0; v < V; v++) {
		// Initialize the distance to all vertices to infinity
		final_dist[v] = INT_MAX;
		// And having a null predecessor
		final_predecessors[v] = -1;
	}
	// The distance from the source to itself is, of course, zero
	final_dist[sourceVertex] = 0;

	// initialize the per thread arrays
	for(int t = 0; t < num_threads; t++) {
		dist_per_thread[t] = (int*) malloc(V * sizeof(int));
		predecessors_per_thread[t] = (int*) malloc(V * sizeof(int));
	}


	// Step 2: relax edges repeatedly
	for (int v = 0; v < V - 1; v++) {

		#pragma omp parallel default(none) \
		shared(V, E, edges, num_threads, dynamic_V, dynamic_E) \
		shared(final_dist, final_predecessors, dist_per_thread, predecessors_per_thread) \
		private(src, dest, weight, dist_src, dist_dest) 
		{
			int id_thread = omp_get_thread_num();

			// update the arrays of each thread
			memcpy(dist_per_thread[id_thread], final_dist, V*sizeof(int));
			memcpy(predecessors_per_thread[id_thread], final_predecessors, V*sizeof(int));
		
			#pragma omp barrier

			#pragma omp for schedule(dynamic, dynamic_E)
			for (int e = 0; e < E; e++) {
				src = edges[e]->src;
				dist_src = dist_per_thread[id_thread][src];

				if (dist_src != INT_MAX){
					dest = edges[e]->dest;
					weight = edges[e]->weight;
					dist_dest = dist_per_thread[id_thread][dest];

					if (dist_src + weight < dist_dest) {
						dist_per_thread[id_thread][dest] = dist_src + weight;
						predecessors_per_thread[id_thread][dest] = src;
					}
				}
			}
		
			#pragma omp barrier
			
			// aggregate the results
			#pragma omp for schedule(dynamic, dynamic_V)
			for (int v1=0; v1 < V; v1++) {
				
				for (int t=0; t < num_threads; t++) {
					if (dist_per_thread[t][v1] < final_dist[v1]){
						final_dist[v1] = dist_per_thread[t][v1];
						final_predecessors[v1] = predecessors_per_thread[t][v1];
					}
						
				}
			}
		}
	}
	

	// Step 3: check for negative-weight cycles
	int negative_cycle = 0;
	#pragma omp parallel for reduction(max:negative_cycle) private(src) private(dest) private(weight) schedule(dynamic, dynamic_E)
    for (int j = 0; j < E; j++) {
		Edge* edge = graph->edges[j];
		src = edge->src;
		dest = edge->dest;
		weight = edge->weight;
        
		if (final_dist[src] != INT_MAX && final_dist[src] + weight < final_dist[dest]) {
			negative_cycle = 1;  
		}
		else{
			negative_cycle = 0;
		}
    }

	if (negative_cycle > 0){
		printf("GRAPH CONTAINS A NEGATIVE-WEIGHT CYCLE\n");
		return;
	}

	free(final_dist);
	free(final_predecessors);
	free(dist_per_thread);
	free(predecessors_per_thread);

    return;  // No negative-weight cycles
}



void parallelBellmanFordNaive(Graph* graph, int sourceVertex) {
    int V = graph->numVertices;
	int E = graph->numEdges;
	Edge** edges = graph->edges;

	Edge* edge;
	int src;
	int dest;
	int weight;

	int num_threads = omp_get_max_threads();
	
	int* dist;
	int* predecessors;
	dist = (int*) malloc(sizeof(int) * V);
	predecessors = (int*) malloc(sizeof(int) * V);

	int dynamic_V = 4 * sqrt(V/num_threads);
	int dynamic_E = 4 * sqrt(E/num_threads);
	dynamic_V = max(1, dynamic_V);
	dynamic_E = max(1, dynamic_E);


	// Step 1: initialize graph
	#pragma omp for schedule(dynamic, dynamic_V)
	for (int i = 0; i < V; i++) {
		// Initialize the distance to all vertices to infinity
		dist[i] = INT_MAX;
		// And having a null predecessor
		predecessors[i] = -1;
	}
	// The distance from the source to itself is, of course, zero
	dist[sourceVertex] = 0;  
	

	// Step 2: relax edges repeatedly
	for (int v = 0; v < V - 1; v++) {

		#pragma omp parallel default(none) \
		shared(V, E, edges, num_threads, dynamic_V, dynamic_E) \
		shared(dist, predecessors) \
		private(src, dest, weight) 
		{
			int id_thread = omp_get_thread_num();
			#pragma omp for schedule(dynamic, dynamic_E)
			for (int e = 0; e < E; e++) {
				dest = edges[e]->dest;
				src = edges[e]->src;
				weight = edges[e]->weight;

				if (dist[src] != INT_MAX) {
					#pragma omp critical
					{
						if (dist[src] + weight < dist[dest]){
							dist[dest] = dist[src] + weight;
							predecessors[dest] = src;
						}
					}
				}
			}
		}
	}

	// Step 3: check for negative-weight cycles
	int negative_cycle = 0;
	#pragma omp parallel for reduction(max:negative_cycle) private(src) private(dest) private(weight) schedule(dynamic, dynamic_E)
    for (int j = 0; j < E; j++) {
		Edge* edge = graph->edges[j];
		src = edge->src;
		dest = edge->dest;
		weight = edge->weight;
        
		if (dist[src] != INT_MAX && dist[src] + weight < dist[dest]) {
			negative_cycle = 1;  
		}
		else{
			negative_cycle = 0;
		}
    }

	if (negative_cycle > 0){
		printf("GRAPH CONTAINS A NEGATIVE-WEIGHT CYCLE\n");
		return;
	}

	free(dist);
	free(predecessors);

    return;  // No negative-weight cycles
}



int main(int argc, char *argv[]){
	// Check for the required argument
	if (argc != 2) {
		printf("Usage: %s <input_file>\n", argv[0]);
		return 1;
	}

	Graph* graph = readFile(argv[1]);

	double tstart, tstop;

	tstart = omp_get_wtime();
	parallelBellmanFordNaive(graph, 0);
	tstop = omp_get_wtime();
	printf("Naive Elapsed time: %f\n", tstop - tstart);
	
	tstart = omp_get_wtime();
	parallelBellmanFord(graph, 0);
	tstop = omp_get_wtime();
	printf("Parallel Elapsed time: %f\n", tstop - tstart);

	tstart = omp_get_wtime();
	parallelBellmanFordMemcpy(graph, 0, 2);
	tstop = omp_get_wtime();
	printf("Memcpy Dynamic Elapsed time: %f\n", tstop - tstart);

	tstart = omp_get_wtime();
	parallelBellmanFordMemcpy(graph, 0, 1);
	tstop = omp_get_wtime();
	printf("Memcpy notDynamic Elapsed time: %f\n", tstop - tstart);

	return 0;
}

