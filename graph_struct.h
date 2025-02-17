#ifndef GRAPH_STRUCT
#define GRAPH_STRUCT

#define MAX_WEIGHT 100  // Maximum weight for an edge
#define MAX_LINE_LENGTH 1024  

#include <stdio.h>
#include <stdlib.h>

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
Graph* createGraph(int numVertices, int numEdges);

// function to add and edge to the graph
void addEdge(Graph* graph, int numEdge, int src, int dest, int weight);

// function to print a graph to a file
void printGraph(Graph* graph, FILE *output_file);

/*
given a file, read a number until there is a \n
it is used to read the number of verteces and edges
*/
int readNumber(FILE *file);

// Given a file containing info, return the graph structure
Graph* readFile(char *fileName);
	
#endif 