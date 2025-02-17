#!/bin/bash
#SBATCH --job-name=generate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mail.com
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=output_generation.txt

# Exit on non-zero return code
set -e

# Compile programs
gcc -fopenmp generate_graph.c -o generate_graph.o -c
gcc -fopenmp -c -o graph_struct.o graph_struct.c
gcc -fopenmp -o generate_graph graph_struct.o generate_graph.o

graph_sizes=(100 1000 2000 5000 10000 20000)

for size in "${graph_sizes[@]}"; do
    ./generate_graph $size "graph$size.txt"
done

