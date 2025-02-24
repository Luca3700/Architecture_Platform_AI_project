# Architecture and Platform for Artificial Intelligence
# Parallel implementation of the Bellman-Ford algorithm

## Introduction
The Bellman-Ford algorithm is an important algorithm that computes shortest paths from a single source vertex to all the others in a weighted directed graph. With respect to other algorithms designed to handle the same problem, like the Dijkstra’s algorithm, it is slower but more versatile, since it is able to handle graphs with one or more negative edge weights. Moreover it will report if there is a negative weight cycle. The presence of the latter indicates that there is no cheapest path to reach a destination vertex: any path that has a point on the negative cycle can be made cheaper by one more walk around the negative cycle.

The Bellman-Ford algorithm has a variety of applications, for example
* network routing: it is used to find the shortest paths in routing tables, helping data packets navigate efficiently across networks
* GPS navigation and logistics: to compute fastest routes between locations, minimizing the cost
* path planning.

In this document parallel implementation of the Bellman-Ford algorithm using two different frame-works are reported:
* OpenMP, it is an application programming interface (API) that supports shared-memory multi-processing programming in different languages and on many platforms
* CUDA, it is a proprietary parallel computing platform and API that allows software to use certain types of graphics processing units (GPUs) for accelerated general-purpose processing

The developed parallel implementations will be analyzed in detail, along with their performance.

## Bellman-Ford Algorithm
