#include "graph_struct.h"

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
