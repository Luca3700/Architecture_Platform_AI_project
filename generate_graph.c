#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "graph_struct.h"

/***
gcc generate_graph_easy.c -o generate_graph_easy.o -c
gcc graph_struct.c -o graph_struct.o -c
gcc -o generate_graph graph_struct.o generate_graph_easy.o
./generate_graph 10 "test.txt"
***/

void addNewEdgeToGraph(Graph* graph, int numEdge){
	int numVertices = graph->numVertices;

	bool isSrcAndDestOk = false;
	int src = 0;
	int dest = 0;

	while (!isSrcAndDestOk){
		src = rand() % numVertices;
    	dest = rand() % numVertices;
    
		// Verify that src and dest are different
		isSrcAndDestOk = (src != dest);

		// Verify that this edge does not exist
		if (isSrcAndDestOk) {
			for (int i=0; i<numEdge; i++){
				int edgeSrc = graph->edges[i]->src;
				int edgeDest = graph->edges[i]->dest;

				if (edgeSrc == src && edgeDest == dest){
					isSrcAndDestOk = false;
				}
			}
		}
	}	

	int weight = (rand() % MAX_WEIGHT) + 1;
	addEdge(graph, numEdge, src, dest, weight);
	// printf("Added: src %d, dest %d, weight %d\n", src, dest, weight);
}


// Helper function for Depth-First Search (DFS) traversal
void DFSUtil(Graph* graph, int vertex, bool visited[]) {
	visited[vertex] = true; // Mark the current vertex as visited

	// Recur for all adjacent vertices
	for (int i = 0; i < graph->numEdges; i++) {
		Edge* edge = graph->edges[i];
		if (edge->src == vertex && !visited[edge->dest]) {
			DFSUtil(graph, edge->dest, visited);
		}
	}
}

// Function to check if all vertices are reachable from the vertex 0
bool isConnected(Graph* graph) {
	int numVertices = graph->numVertices;
	bool visited[numVertices]; // Array to keep track of visited vertices

	// Initialize all visited to false
	for (int i = 0; i < numVertices; i++) {
		visited[i] = false;
	}

	// Call a helper function to perform DFS starting from vertex 0
	DFSUtil(graph, 0, visited);

	// Check if all vertices are marked as visited
	for (int i = 0; i < numVertices; i++) {
		if (!visited[i]) {
			return false; // If any vertex is not visited, graph is not connected
		}
	}
	return true; // If all vertices are visited, graph is connected
}


void generateRandomGraph(Graph* graph){
  srand(time(NULL)); // set the random seed

  int numVertices = graph->numVertices;
  int numEdges = graph->numEdges;

  int count = 0;
  bool connected = true;

  while(count < numEdges){ // || !connected){
	addNewEdgeToGraph(graph, count);
    count++;

    // if (count >= numEdges){
    //   //Checking if it is connected
    //   connected = isConnected(graph);
    //   // printf("\nConnected? %d\n", connected);
    //   if (!connected){
    //     // Increase the size of the edges array
    //     Edge** new_edges = (Edge**)realloc(graph->edges, (graph->numEdges + 1) * sizeof(Edge*));
    //     if (new_edges == NULL) {
    //       printf("Memory allocation failed!\n");
    //       exit(1);
    //     }
    //     graph->edges = new_edges;
     
    //     // Increment the edge count
    //     graph->numEdges++;
    //   }
    // }
  }
}



int main(int argc, char *argv[]){
	// Check for the required argument
	if (argc != 3) {
		printf("Usage: %s <numVertices> <output_file>\n", argv[0]);
	return 1;
	}

	// Open the output file
	FILE *output_file = fopen(argv[2], "w");
	if (output_file == NULL) {
		printf("Error opening output file!\n");
		return 1;
	}

	// Read the numVertices from the command-line argument
	int numVertices = atoi(argv[1]);
	// int numEdges = numVertices/2;
	int numEdges = numVertices * 10;
 	
	// printf("Creating the graph\n");
	Graph* graph = createGraph(numVertices, numEdges);

	generateRandomGraph(graph);
	
	// Redirect output to the file
	// printf("Generated graph:\n");
 	printGraph(graph, output_file);  // Output will be written to the file

	fclose(output_file);  

	return 0;
}

