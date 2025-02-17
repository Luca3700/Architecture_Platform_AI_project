#!/bin/bash
#SBATCH --job-name=graph1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luca.torzi@studio.unibo.it
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=output.txt
#SBATCH --gres=gpu:1

# Exit on non-zero return code
set -e

# Compile programs
gcc -fopenmp generate_graph.c -o generate_graph.o -c
gcc -fopenmp -c -o graph_struct.o graph_struct.c
gcc -fopenmp -o generate_graph graph_struct.o generate_graph.o

gcc -fopenmp -c -o BellmanFord.o BellmanFord.c
gcc -o BellmanFord graph_struct.o BellmanFord.o -fopenmp -lm

nvcc BellmanFord.cu -o BellmanFordCuda
nvcc BellmanFord_blkdim1.cu -o BellmanFordCuda_blkdim1

# Array of graph sizes to test
graph_sizes=(100 1000 2000 5000 10000 20000)

# Function to run BellmanFord with different thread counts
run_bellmanford() {
    local threads=$1
    local size=$2
    echo "BellmanFord $threads"
    for i in {1..10}; do
        OMP_NUM_THREADS=$threads ./BellmanFord "graph$size.txt"
    done
}

# Function to run CUDA versions
run_cuda() {
    local size=$1
    echo "Cuda"
    for i in {1..10}; do
        ./BellmanFordCuda "graph$size.txt"
    done
    echo "Cuda blkdim 1"
    for i in {1..10}; do
        ./BellmanFordCuda_blkdim1 "graph$size.txt"
    done
}

# Main loop for each graph size
for size in "${graph_sizes[@]}"; do
    #./generate_graph $size "graph.txt"
    echo "numVertices: $size"
    echo -e "numEdges: $(expr $size \* 10)\n"

    # Run OpenMP versions with different thread counts
    for threads in 1 2 4 8 16 32 64; do
        run_bellmanford $threads $size
    done

    # Run CUDA versions
    run_cuda $size

    echo "----------------------------------------"
done

# Clean up
rm generate_graph.o graph_struct.o generate_graph BellmanFord.o BellmanFord BellmanFordCuda BellmanFordCuda_blkdim1

