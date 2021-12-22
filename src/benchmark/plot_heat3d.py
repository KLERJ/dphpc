import os
import matplotlib.pyplot as plt
import numpy as np

# Set these for your machine:
rootdir = '/Users/lukevolpatti/Downloads/heat3d'
outdir  = '/Users/lukevolpatti/Downloads/plots'

# Parameters
programs = ['mpi', 'mpi_avx2', 'omp']

def plot_by_size(subfolder):
    subfolder_info = subfolder.split("_")
    num_iters = subfolder_info[-1][:-1]
    prob_size = subfolder_info[-2]

    mpi_total_runtime = {}
    mpi_avx2_total_runtime = {}
    omp_total_runtime = {}
    mpi_per_iter = {}
    avx2_per_iter = {}
    omp_per_iter = {}
    instance_dirs = [ f.path for f in os.scandir(subfolder) if f.is_dir() ]

    for instance_dir in instance_dirs:
        instance_dir_info = instance_dir.split('_')
        nprocs = int(instance_dir_info[-2][:-1])
        program_type = instance_dir_info[-3]

        print(program_type, nprocs)
        
        # Get total runtime from output.txt
        output_txt = os.path.join(instance_dir, 'output.txt')
        with open(output_txt, 'r') as f:
            lines = f.read().splitlines()
            total_runtime = 0

            # TODO: fix lame heuristic for figuring out if there's data or not
            if len(lines) == 0 or len(lines[-1]) > 15:
                total_runtime = 0
            else:
                total_runtime = float(lines[-1])

            if program_type == "mpi":
                mpi_total_runtime[nprocs] = float(total_runtime)
            elif program_type == "avx2":
                mpi_avx2_total_runtime[nprocs] = float(total_runtime)
            elif program_type == "omp":
                omp_total_runtime[nprocs] = float(total_runtime)

        # Get per processor per iteration runtimes
        benchmark_files = [ f.path for f in os.scandir(instance_dir) if f.is_file() ]
        for benchmark_file in benchmark_files:
            if 'benchmark_' in benchmark_file:
                with open(benchmark_file, 'r') as f:
                    values = f.read().splitlines()[-1].split(',')
                    values = np.array([float(i) for i in values])

                    if program_type == "mpi":
                        if nprocs in mpi_per_iter:
                            mpi_per_iter[nprocs] = np.append(mpi_per_iter[nprocs], values)
                        else:
                            mpi_per_iter[nprocs] = values
            
    print(mpi_total_runtime)
    
    # Total runtime plots
    mpi_l = sorted(mpi_total_runtime.items())
    mpi_x, mpi_y = zip(*mpi_l)
    avx2_l = sorted(mpi_avx2_total_runtime.items())
    avx2_x, avx2_y = zip(*avx2_l)
    omp_l = sorted(omp_total_runtime.items())
    omp_x, omp_y = zip(*omp_l)
    plt.figure()
    plt.plot(mpi_x, mpi_y, label='mpi')
    plt.plot(avx2_x, avx2_y, label='avx2')
    plt.plot(omp_x, omp_y, label='omp')
    plt.legend()
    plt.xlabel("# cores")
    plt.ylabel("Total runtime (s)")
    plt.title("Problem size: " + prob_size + ", iterations: " + num_iters)
    plt.show()
    
    # Per processor per iteration plots
    mpi_l = sorted(mpi_per_iter.items())
    mpi_x, mpi_y = zip(*mpi_l)
    mpi_err = [np.std(i) for i in mpi_y]
    mpi_y = [np.mean(i) for i in mpi_y]
    plt.figure()
    plt.errorbar(mpi_x, mpi_y, yerr=mpi_err, ecolor='black', label='mpi')
    plt.legend()
    plt.xlabel("# cores")
    plt.ylabel("Time per iteration per core (s)")
    plt.title("Problem size: " + prob_size + ", iterations: " + num_iters)
    plt.show()

    
if __name__ == '__main__':
    subfolders = [ f.path for f in os.scandir(rootdir) if f.is_dir() ]
    print(subfolders)

    for subfolder in subfolders:
        plot_by_size(subfolder)

