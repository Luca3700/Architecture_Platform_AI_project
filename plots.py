import json

device1 = "cpu"
device2 = "gpu"

def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()

data = read_file("output.txt")

res = {}

sections = data.split("----------------------------------------")

for section in sections:
    lines = section.split("\n")
    numVertices = -1
    numEdges = -1

    device = " "
    num_thread = -1
    elapsed_time = 0
    elapsed_time_dict = {"n":0, "p":0, "m":0, "mn":0}

    current_check = " "

    if len(lines) > 10:

        for line in lines:

            if line.startswith("numVertices"):
                numVertices = int(line.split(':')[1].strip())
                res[numVertices] = {device1:{}, device2:{}}
            elif line.startswith("numEdges"):
                numEdges = int(line.split(':')[1].strip())
            elif line.startswith("BellmanFord") or line == "Cuda" or line == "Cuda blkdim 1":
                # save the old results
                if current_check == "B":
                    elapsed_time_dict["n"] = elapsed_time_dict["n"] / 10
                    elapsed_time_dict["p"] = elapsed_time_dict["p"] / 10
                    elapsed_time_dict["m"] = elapsed_time_dict["m"] / 10
                    elapsed_time_dict["mn"] = elapsed_time_dict["mn"] / 10
                    res[numVertices][device][num_thread] = elapsed_time_dict
                elif current_check == "C":
                    res[numVertices][device][num_thread] = elapsed_time/10

                elapsed_time = 0
                elapsed_time_dict = {"n":0, "p":0, "m":0, "mn":0}

                if line.startswith("BellmanFord"):
                    device = device1
                    num_thread = line.split(" ")[1]
                    current_check = "B"
                elif line == "Cuda" or line == "Cuda blkdim 1":
                    device = device2
                    if len(line.split(" ")) > 1:
                        num_thread = 1
                    else:
                        num_thread = 256
                    current_check = "C"

            elif 'Elapsed time' in line:

                if "Cuda Elapsed time" in line or "Cuda 1 blkdim Elapsed time" in line:
                    seconds = float(line.split(":")[1].strip().split(" ")[0])
                    #print(seconds)
                    elapsed_time += seconds
                    # print(elapsed_time)

                else:
                    seconds = float(line.split(":")[1].strip())
                    if "Naive" in line:
                        elapsed_time_dict["n"] = elapsed_time_dict["n"] + seconds
                    elif "Parallel" in line:
                        elapsed_time_dict["p"] = elapsed_time_dict["p"] + seconds
                    elif "Memcpy Dynamic" in line: 
                        elapsed_time_dict["m"] = elapsed_time_dict["m"] + seconds
                    elif "Memcpy notDynamic" in line: 
                        elapsed_time_dict["mn"] = elapsed_time_dict["mn"] + seconds

        # saving the last result
        res[numVertices][device][num_thread] = elapsed_time/10

# print(res)

with open('output.json', 'w') as f:
    json.dump(res, f, indent=4)



import json
import matplotlib.pyplot as plt

# Load the data
with open('output.json', 'r') as f:
    data = json.load(f)

# Define the number of threads and vertex counts
threads = [1, 2, 4, 8, 16, 32, 64]
vertex_counts = ['100', '1000', '2000', '5000', '10000', '20000']
execution_types = ['p', 'm', 'n', 'mn']
execution_types_translation = {'n':'naive', 'p':'parallel', 'm':'memcpy', 'mn':'memcpy, not dynamic chunk size'}

# Set font size for all plots
plt.rcParams.update({'font.size': 18})

# Function to compute speedup
def compute_speedup(data, vertex_count, execution_type):
    base_time = data[vertex_count]['cpu']['1'][execution_type]
    return [base_time / data[vertex_count]['cpu'][str(t)][execution_type] for t in threads]

# Function to compute strong scaling efficiency
def compute_efficiency(speedup):
    return [s / t for s, t in zip(speedup, threads)]

# Create speedup plots
plt.figure(figsize=(20, 15))
for i, exec_type in enumerate(execution_types):
    plt.subplot(2, 2, i+1)
    for vertex_count in vertex_counts:
        speedup = compute_speedup(data, vertex_count, exec_type)
        plt.plot(threads, speedup, marker='o', label=f'{vertex_count} vertices')
    plt.xscale('log', base=2)
    plt.xlabel('Number of Threads')
    plt.ylabel('Speedup')
    plt.title(f'Speedup Analysis - {execution_types_translation[exec_type]}')
    plt.legend()
    plt.grid(True)
plt.tight_layout()
plt.savefig('speedup_analysis.png')
plt.close()

# Create strong scaling efficiency plots
plt.figure(figsize=(20, 8))
for i, exec_type in enumerate(execution_types[0:2]):
    plt.subplot(1, 2, i+1)
    for vertex_count in vertex_counts:
        speedup = compute_speedup(data, vertex_count, exec_type)
        efficiency = compute_efficiency(speedup)
        plt.plot(threads, efficiency, marker='o', label=f'{vertex_count} vertices')
    plt.xscale('log', base=2)
    plt.xlabel('Number of Threads')
    plt.ylabel('Strong Scaling Efficiency')
    plt.title(f'Strong Scaling Efficiency - {execution_types_translation[exec_type]}')
    plt.legend()
    plt.grid(True)
plt.tight_layout()
plt.savefig('strong_scaling_efficiency.png')
plt.close()

# Create execution time plot
plt.figure(figsize=(20, 8))
min_time = float('inf')
max_time = 0
for i, exec_type in enumerate(execution_types[0:2]):
    plt.subplot(1, 2, i+1)
    for vertex_count in vertex_counts:
        times = [data[vertex_count]['cpu'][str(t)][exec_type] for t in threads]
        min_time = min(min_time, min(times))
        max_time = max(max_time, max(times))
        plt.plot(threads, times, marker='o', label=f'{vertex_count} vertices')
    plt.xscale('log', base=2)
    plt.yscale('log', base=10)
    plt.xlabel('Number of Threads')
    plt.ylabel('Elapsed Time (s)')
    plt.title(f'Elapsed Time - CPU {execution_types_translation[exec_type]}')
    plt.legend()
    plt.grid(True)

# Set the same y-axis limit for all subplots
for ax in plt.gcf().get_axes():
    ax.set_ylim(min_time * 0.9, max_time * 1.1)

plt.tight_layout()
plt.savefig('elapsed_time_cpu.png')
plt.close()

# Create elapsed time plot for GPU
plt.figure(figsize=(10, 7))
for vertex_count in vertex_counts:
    times = [data[vertex_count]['gpu']['1'], data[vertex_count]['gpu']['256']]
    plt.plot(['1', '256'], times, marker='o', label=f'{vertex_count} vertices')
plt.yscale('log', base=10)
plt.xlabel('Number of Threads')
plt.ylabel('Elapsed Time (s)')
plt.title('Elapsed Time - GPU')
plt.legend()
# plt.grid(True)
plt.tight_layout()
plt.savefig('elapsed_time_gpu.png')
plt.close()

# Create plot for elapsed time vs number of vertices on a logarithmic scale
vertex_count_int = [int(vc) for vc in vertex_counts]  # Convert vertex counts to integers for plotting

cpu_times = [data[vc]['cpu']['4']['m'] for vc in vertex_counts]
gpu_times1 = [data[vc]['gpu']['1'] for vc in vertex_counts]
gpu_times2 = [data[vc]['gpu']['256'] for vc in vertex_counts]

plt.figure(figsize=(12, 8))
plt.plot(vertex_count_int, cpu_times, marker='o', label='CPU (4 threads, memcpy version)')
plt.plot(vertex_count_int, gpu_times1, marker='o', label='GPU (1 thread per block)')
plt.plot(vertex_count_int, gpu_times2, marker='o', label='GPU (256 threads per block)')

plt.xscale('log', base=10)
plt.yscale('log', base=10)
plt.xlabel('Number of Vertices')
plt.ylabel('Elapsed Time (s)')
plt.title('Elapsed Time of the CUDA implementation w.r.t. OpenMP')
plt.legend()
# plt.grid(True)
plt.tight_layout()
plt.savefig('elapsed_time_vs_vertices.png')
plt.close()

