#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#include <cuda.h> 
#include <cuda_runtime_api.h> 
#include <cuda_runtime.h>


#define MAX_WEIGHT 100  // Maximum weight for an edge
#define MAX_LINE_LENGTH 1024  

typedef struct {
  int src; // source vertex of the edge
  int dest;  // destination vertex of the edge
  int weight;  // weight of the edge
} Edge;

typedef struct {
  int numVertices;
  int numEdges;
  Edge** edges;  // array of edges
} Graph;

// function to create an empty graph
Graph* createGraph(int numVertices, int numEdges) {
	//creating the graph
	Graph* graph = (Graph*)malloc(sizeof(Graph));
	graph->numVertices = numVertices;
	graph->numEdges = numEdges;

	// Creating a list of edges
	graph->edges = (Edge**)malloc(numEdges * sizeof(Edge));

	return graph;
}

// function to add and edge to the graph
void addEdge(Graph* graph, int numEdge, int src, int dest, int weight){
	// Allocate the new edge
	graph->edges[numEdge] = (Edge*)malloc(sizeof(Edge));

	Edge* newEdge = graph->edges[numEdge];
	// insert the information of the edge
	newEdge->src = src;
	newEdge->dest = dest;
	newEdge->weight = weight;
}

// function to print a graph to a file
void printGraph(Graph* graph, FILE *output_file){
	int numVertices = graph->numVertices;
	int numEdges = graph->numEdges;

	fprintf(output_file, "%d\n", numVertices);
	fprintf(output_file, "%d\n", numEdges);
	
	printf("numVertices: %d\n", numVertices);
	printf("numEdges: %d\n\n", numEdges);

	// for each vertex
	for (int v=0; v<numVertices; v++) {
		fprintf(output_file, "%d:", v);
		// printf("%d:", v);

		// for each edges
		for (int e=0; e<numEdges; e++){
			Edge* edge = graph->edges[e];
			if (edge->src == v){
				fprintf(output_file, "%d ", edge->dest);
				// printf("%d ", edge->dest);
				fprintf(output_file, "%d,", edge->weight);
				// printf("%d,", edge->weight);
			}
		}
		fprintf(output_file, "\n");
		//printf("\n");  
	}
}

int readNumber(FILE *file){
	/*
	given a file, read a number until there is a \n
	it is used to read the number of verteces and edges
	*/
	char c;
	char *string;
	int n = 0;
	int res;

	string = (char*) malloc(sizeof(char) * MAX_LINE_LENGTH);
	
	// read the character until end of line
	c = (char) fgetc(file);
	while (c != '\n'){
		string[n] = c;
		n++;
		c = (char) fgetc(file);
	}
	string[n] = '\0';
	res = atoi(string);
	
	// free the memory and return
	free(string);
	return res;
}

Graph* readFile(char *fileName){
	/*
	Given a file containing info, return the graph structure
	*/
	FILE *file = fopen(fileName, "r");
	char *string;
	int n = 0;
    
	int c1; // character (int)
	char c2; // character (char)

	int numVertices;
	int numEdges;

	int srcVertex;
	int destVertex;
	int weightEdge;

	int count = 0;

    if (file == NULL){
        return NULL; //could not open file
	}

	// printf("Reading the number of vertices\n");
	numVertices = readNumber(file);
	// printf("Reading the number of edges\n");
	numEdges = readNumber(file);

	// creating an empty graph
	// printf("Creating an empty graph\n");
	Graph* graph = createGraph(numVertices, numEdges);

	string = (char*) malloc(sizeof(char) * MAX_LINE_LENGTH);

	c1 = fgetc(file);
	while (c1 != EOF){
		c2 = (char) c1;

		//read the src vertex	
		while (c2 != ':'){
			string[n++] = c2;
			c1 = fgetc(file);
			c2 = (char) c1;
		}
		string[n] = '\0';
		srcVertex = atoi(string);
			
		n = 0;
		// read all the destinations
		c1 = fgetc(file);
		c2 = (char) c1;

		while (c2 != '\n'){

			// read the destVertex
			while (c2 != ' '){
				string[n++] = c2;
				c1 = fgetc(file);
				c2 = (char) c1;
			}
			string[n] = '\0';
			destVertex = atoi(string);

			// read the weight of the edge
			n = 0;
			c1 = fgetc(file);
			c2 = (char) c1;
			while (c2 != ','){
				string[n++] = c2;
				c1 = fgetc(file);
				c2 = (char) c1;
			}
			string[n] = '\0';
			weightEdge = atoi(string);

			// add the edge to the graph
			addEdge(graph, count++, srcVertex, destVertex, weightEdge);

			// restarting the loop to search other edges with the same source
			n = 0;
			c1 = fgetc(file);
			c2 = (char) c1;
		}
		c1 = fgetc(file);
    }

	free(string);  
	fclose(file);   

	return graph;
}

void printArray(int* array, int n){
	for (int i=0; i<n; i++){
		printf("%d: %d\n", i, array[i]);
	}
}













#define BLKDIM 256

__device__ long long pack(int dist, int pred) {
    return ((long long)dist << 32) | (pred & 0xFFFFFFFF);
}

__device__ __host__ int unpackDist(long long val) {
    return (int)(val >> 32);
}

__device__ __host__ int unpackPred(long long val) {
    return (int)(val & 0xFFFFFFFF);
}


__global__ void step1_init(long long* dist_pred, int numVertices, int sourceVertex){
	const int idx = blockIdx.x * BLKDIM + threadIdx.x;

	if (idx < numVertices) {
		// Initialize the distance to all vertices to infinity and predecessor to -1
		dist_pred[idx] = pack(INT_MAX, -1);
	}

	// The distance from the source to itself is zero
	if (idx == sourceVertex){
		dist_pred[idx] = pack(0, -1);
	}
}

__global__ void relaxEdges(int* srcs, int* dests, int* weights, long long* dist_pred, int numEdges){
	const int tid = threadIdx.x;
	const int idx = blockIdx.x * BLKDIM + tid;

	if (idx < numEdges){
		int src = srcs[idx];
		// Extract distance from the source
		int dist_src = unpackDist(dist_pred[src]);

		if (dist_src != INT_MAX) {
			int dest = dests[idx];
			int weight = weights[idx];
			int new_dist = dist_src + weight;

			// Create the long long variable
			long long new_dist_pred = pack(new_dist, src);

			if (new_dist < unpackDist(dist_pred[dest])) {
				// Update the array
				atomicMin(&dist_pred[dest], new_dist_pred);
			}
		}
	}
}

__global__ void check_negative_cycle(int* srcs, int* dests, int* weights, long long* dist_pred, int numEdges, int* negative_cycles){
	__shared__ int local_res[BLKDIM];
	const int idx = blockIdx.x * BLKDIM + threadIdx.x;
	const int tid = threadIdx.x;

	if (idx < numEdges){ 
		int src = srcs[idx];
		int dest = dests[idx];
		int weight = weights[idx];
		int dist_src = unpackDist(dist_pred[src]);

		if (dist_src != INT_MAX && dist_src + weight < unpackDist(dist_pred[dest])) {
			local_res[tid] = 1;  
		}
		else{
			local_res[tid] = 0;
		}
	}
	else{
		local_res[tid] = 0;
	}

	__syncthreads();

	// reduction inside the block
    for (int stride = BLKDIM / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            local_res[tid] += local_res[tid + stride];
        }
        __syncthreads();
    }

	// saving the final value
	negative_cycles[blockIdx.x] = local_res[0];
}

void cudaParallelBellmanFord(Graph* graph, int* srcs, int* dests, int* weights, int sourceVertex){
	int V = graph->numVertices;
	int E = graph->numEdges;

	const size_t size_V = V * sizeof(long long);
	const size_t size_E = E * sizeof(int);

	// initializing the results of the algorithm
	long long* dist_pred; 
	long long* d_dist_pred;
	dist_pred = (long long*) malloc(size_V);
	cudaMalloc((void **)&d_dist_pred, size_V);

	int how_many_blocks_V = (V + BLKDIM - 1) / BLKDIM;
	int how_many_blocks_E = (E + BLKDIM - 1) / BLKDIM;

	// initializing the copies of the edge information on the device	
	int* d_srcs; 
	int* d_dests;
	int* d_weights;
	cudaMalloc((void **)&d_srcs, size_E);
	cudaMalloc((void **)&d_dests, size_E);
	cudaMalloc((void **)&d_weights, size_E);
	cudaMemcpy(d_srcs, srcs, size_E, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dests, dests, size_E, cudaMemcpyHostToDevice);
	cudaMemcpy(d_weights, weights, size_E, cudaMemcpyHostToDevice);
	cudaError_t status = cudaGetLastError();
	if (status != cudaSuccess) {
		printf("Error: %s launching kernel\n", cudaGetErrorString(status));
		exit(-1);
	} 

	// Step 1: initialize graph
	step1_init<<<how_many_blocks_V, BLKDIM>>>(d_dist_pred, V, sourceVertex);
	status = cudaGetLastError();
	if (status != cudaSuccess) {
		printf("Error: %s launching kernel\n", cudaGetErrorString(status));
		exit(-1);
	} 

	// Step 2: relax edges repeatedly
	for (int v = 0; v < V - 1; v++) {
		relaxEdges <<<how_many_blocks_E, BLKDIM>>> (d_srcs, d_dests, d_weights, d_dist_pred, E);
		cudaDeviceSynchronize();
		
		status = cudaGetLastError();
		if (status != cudaSuccess) {
			printf("Error: %s launching kernel\n", cudaGetErrorString(status));
			exit(-1);
		} 
	}
	
	// Step 3: check for negative-weight cycles
	int* negative_cycles;
	int* d_negative_cycles;
	negative_cycles = (int*)malloc(how_many_blocks_E * sizeof(int));
	cudaMalloc((void **)&d_negative_cycles, how_many_blocks_E * sizeof(int));

	check_negative_cycle <<<how_many_blocks_E, BLKDIM>>> (d_srcs, d_dests, d_weights, d_dist_pred, E, d_negative_cycles);
	status = cudaGetLastError();
	if (status != cudaSuccess) {
		printf("Error: %s launching kernel\n", cudaGetErrorString(status));
		exit(-1);
	} 
	// Copy back the array containing the result of the computation
	cudaMemcpy(negative_cycles, d_negative_cycles, how_many_blocks_E * sizeof(int), cudaMemcpyDeviceToHost);

	// reduction of the values returned by each block
	int negative_cycle = 0;
	for (int num_block = 0; num_block < how_many_blocks_E; num_block++) {
		negative_cycle += negative_cycles[num_block];
	}
	if (negative_cycle > 0) {
		printf("GRAPH CONTAINS A NEGATIVE-WEIGHT CYCLE\n");
	}

	// copy back the results
	cudaMemcpy(dist_pred, d_dist_pred, size_V, cudaMemcpyDeviceToHost);

	int* dist = (int*)malloc(V * sizeof(int));
	int* predecessors = (int*)malloc(V * sizeof(int));
	for (int i = 0; i < V; i++) {
		dist[i] = unpackDist(dist_pred[i]);
		predecessors[i] = unpackPred(dist_pred[i]);
	}

	cudaFree(d_dist_pred);
	cudaFree(d_srcs);
	cudaFree(d_dests);
	cudaFree(d_weights);
	cudaFree(d_negative_cycles);

    // printf("Distances (from: distance):\n");
	// printArray(dist, V);
	// printf("Predecessors (of: predecessor):\n");
	// printArray(predecessors, V);

	free(dist_pred);
	free(negative_cycles);
	free(dist);
	free(predecessors);

	return;
}

int main(int argc, char *argv[]) {
	
	if (argc != 2) {
		printf("Usage: %s <input_file>\n", argv[0]);
		return 1;
	}

	Graph* graph = readFile(argv[1]);
	
	Edge** edges = graph->edges;
	int numEdges = graph->numEdges;

	int* srcs = (int*)malloc(numEdges * sizeof(int));
	int* dests = (int*)malloc(numEdges * sizeof(int));
	int* weights = (int*)malloc(numEdges * sizeof(int));

	for (int i = 0; i < numEdges; i++){
		srcs[i] = edges[i]->src;
		dests[i] = edges[i]->dest;
		weights[i] = edges[i]->weight;
	}

	float elapsed_time = -1;
	cudaEvent_t tstart, tstop;

	cudaEventCreate(&tstart);
	cudaEventCreate(&tstop);
	cudaEventRecord(tstart, 0);

	cudaParallelBellmanFord(graph, srcs, dests, weights, 0);
	
	cudaEventRecord(tstop, 0);
	cudaEventSynchronize(tstop);
	cudaEventElapsedTime(&elapsed_time, tstart, tstop);   
	
	elapsed_time /= 1000;
	printf("Cuda Elapsed time: %f s\n", elapsed_time);

	cudaEventDestroy(tstart);
	cudaEventDestroy(tstop);

	return 0;
}
