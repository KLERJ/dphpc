#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse
import os
import numpy as np

# ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument('results_path', help='root directory of the benchmark results')
parser.add_argument('dim',type=int,  help='input size, 2^dim x 2^dim')
args = parser.parse_args()
dim = args.dim
print('dim=', dim)

########### GRAPH PROPERTIES ############

# handle weak scaling dimensions
dim_w = dim % 100
W_POW_OF_2 = '$2^{' + str(dim_w) +  '}$'
if dim > 100:
    dim_h = dim // 100
else:
    dim_h = dim_w
H_POW_OF_2 = '$2^{' + str(dim_h) +  '}$'
X_LABEL = '# cores'
FILENAME = 'deriche_dim' + str(dim)
AXIS_LABEL_COLOR = (.4, .4, .4)
BACKGROUND_COLOR = '#eeeeee'
HORIZONTAL_LINES_COLOR = '#fafafa'

targets = ['deriche_mpi_baseline', 'deriche_mpi_rdma', 'deriche_mpi_rdma', 'deriche_omp', 'deriche_ref', 'deriche_mpi_segments', 'deriche_mpi_segments']
sws = [1, 1, 16, 1, 1, 32, 1]
labels = ['MPI Alltoall', 'MPI RDMA SW = 1',  'MPI RDMA SW = 16' , 'OpenMP', 'polybench', 'MPI Segments SW = 32', 'MPI Segments SW = 1']


n_cpus = [1, 2, 4, 8, 16, 32]
x_ticks = [x for x in range(6)]
x_labels = [str(2**x) for x in range(6)]
line_styles = ['b.-', 'r.-', 'g.-','c.-',  'y.-', 'm.-', 'k.-']

# dims = [10, 11, 12, 13, 14]

################# DATA ##################

def read_results_default(root_dir_path):
    results = {}
    for dim_dir in os.listdir(root_dir_path):
        dim_dir_path = os.path.join(root_dir_path, dim_dir)
        dim = int(dim_dir.split('_')[1])
        dim_results = {}
        for cpu_dir in os.listdir(dim_dir_path):
            cpu_dir_path = os.path.join(dim_dir_path, cpu_dir)
            cpus = int(cpu_dir.split('_')[1])
            cpu_results = {}
            for rep_dir in os.listdir(cpu_dir_path):
                rep_dir_path = os.path.join(cpu_dir_path, rep_dir)
                rep = int(rep_dir.split('_')[1])
                # print(rep_dir_path," ",n_cpus)
                bm_results = []
                for i in range(cpus):
                    bm_file_path = os.path.join(rep_dir_path, 'benchmark_' + str(i) + '.csv')
                    # Read file into list, append list to bm_results
                    bm_file_res = []
                    file = open(bm_file_path, 'r')
                    for line in file:
                        bm_file_res.append(float(line))
                    bm_results.append(bm_file_res[3])
                # Add results of this rep to results of this cpu number
                cpu_results[rep] = max(bm_results)
            # Add results of this number of cpu to results of this dim
            dim_results[cpus] = cpu_results
        results[dim] = dim_results

    return results

def read_results_segments(root_dir_path):
    results = {}
    for sw_dir in os.listdir(root_dir_path):
        sw_results = {}
        sw_dir_path = os.path.join(root_dir_path, sw_dir)
        sw = int(sw_dir.split('_')[1])
        for dim_dir in os.listdir(sw_dir_path):
            dim_dir_path = os.path.join(sw_dir_path, dim_dir)
            dim = int(dim_dir.split('_')[1])
            dim_results = {}
            for cpu_dir in os.listdir(dim_dir_path):
                cpu_dir_path = os.path.join(dim_dir_path, cpu_dir)
                cpus = int(cpu_dir.split('_')[1])
                cpu_results = {}
                for rep_dir in os.listdir(cpu_dir_path):
                    rep_dir_path = os.path.join(cpu_dir_path, rep_dir)
                    rep = int(rep_dir.split('_')[1])
                    # print(rep_dir_path," ",n_cpus)
                    bm_results = []
                    for i in range(cpus):
                        bm_file_path = os.path.join(rep_dir_path, 'benchmark_' + str(i) + '.csv')
                        # Read file into list, append list to bm_results
                        bm_file_res = []
                        # try:
                        file = open(bm_file_path, 'r')
                        # except FileNotFoundError:
                            # print(f'file {file} not found - ignoring error')
                            # continue
                        for line in file:
                            bm_file_res.append(float(line))
                        bm_results.append((bm_file_res[3]))
                    # Add results of this rep to results of this cpu number
                    cpu_results[rep] = max(bm_results)
                # Add results of this number of cpu to results of this dim
                dim_results[cpus] = cpu_results
            sw_results[dim] = dim_results
        results[sw] = sw_results

    return results

def read_results_omp(root_dir_path):
    results = {}
    for dim_dir in os.listdir(root_dir_path):
        dim_dir_path = os.path.join(root_dir_path, dim_dir)
        dim = int(dim_dir.split('_')[1])
        dim_results = {}
        for cpu_dir in os.listdir(dim_dir_path):
            cpu_dir_path = os.path.join(dim_dir_path, cpu_dir)
            cpus = int(cpu_dir.split('_')[1])
            cpu_results = {}
            for rep_dir in os.listdir(cpu_dir_path):
                rep_dir_path = os.path.join(cpu_dir_path, rep_dir)
                rep = int(rep_dir.split('_')[1])
                bm_results = []
                bm_file_path = os.path.join(rep_dir_path, 'stdout')
                # Read file into list, append list to bm_results
                file = open(bm_file_path, 'r')
                for line in file:
                    bm_results.append(float(line))
                # Add results of this rep to results of this cpu number
                cpu_results[rep] = bm_results
            # Add results of this number of cpu to results of this dim
            dim_results[cpus] = cpu_results
        results[dim] = dim_results
    return results

def read_results_ref(root_dir_path):
    results = {}
    for dim_dir in os.listdir(root_dir_path):
        dim_dir_path = os.path.join(root_dir_path, dim_dir)
        dim = int(dim_dir.split('_')[1])
        dim_results = {}
        for rep_dir in os.listdir(dim_dir_path):
            rep_dir_path = os.path.join(dim_dir_path, rep_dir)
            rep = int(rep_dir.split('_')[1])
            # print(rep_dir_path," ",n_cpus)
            bm_results = []
            bm_file_path = os.path.join(rep_dir_path, 'stdout')
            # Read file into list, append list to bm_results
            file = open(bm_file_path, 'r')
            for line in file:
                bm_results.append(float(line))
            # Add results of this rep to results of this dim
            # print(bm_file_path)
            # print(bm_results)
            dim_results[rep] = bm_results[0]
        results[dim] = dim_results
    return results

def runtime_for_dim(sw, target):
    # print('target: ', target)

    root_dir_path = os.path.join(args.results_path, target)
    # whether the results contains a 'SW' parameter or not
    is_segmented = (target == 'deriche_mpi_segments') | (target == 'deriche_mpi_rdma')
    is_mpi = is_segmented | (target == 'deriche_mpi_baseline')
    is_omp = (target == 'deriche_omp')

    results = {}
    if is_segmented:
        # print('sw: ', sw)
        results = read_results_segments(root_dir_path)
    elif is_mpi:
        results = read_results_default(root_dir_path)
    elif is_omp:
        results = read_results_omp(root_dir_path)
    elif target == 'deriche_ref':
        results = read_results_ref(root_dir_path)[dim]
        dim_baseline_reps = np.array([rep for rep in results.values()])
        mean = [ np.mean(dim_baseline_reps) for n in range(len(n_cpus)) ]
        var = [ np.var(dim_baseline_reps) for n in range(len(n_cpus))]
        return mean, var
    else:
        print("Invalid target parameter")
        exit(-1)

    dim_results = {}
    if is_segmented:
        dim_results = results[sw][dim]
    else:
        dim_results = results[dim]

    # Compute arithmetic mean and variance, speedup compared to baseline, efficiency
    dim_cpus_reps = np.array([[i for i in dim_results[cpu].values()] for cpu in n_cpus])
    mean = np.array([ np.mean(reps) for reps in dim_cpus_reps ])
    variance = np.array([ np.var(reps) for reps in dim_cpus_reps ])

    return mean, variance


# Weak scaling speedup data
def speedup_for_dims(sw, target, dims):

    baseline_root_path = os.path.join(args.results_path, 'deriche_ref')
    results_ref = read_results_ref(baseline_root_path)
    baseline_ws = []
    for d in dims:
        baseline = results_ref[d]
        dim_baseline_reps = np.array([rep for rep in baseline.values()])
        dim_baseline_mean = np.mean(dim_baseline_reps)
        dim_baseline_var = np.var(dim_baseline_reps)
        baseline_ws.append(dim_baseline_mean)

    print('baseline_ws: ', baseline_ws)

    root_dir_path = args.results_path + '/' + target
    # whether the results contains a 'SW' parameter or not
    is_segmented = (target == 'deriche_mpi_segments') | (target == 'deriche_mpi_rdma')
    is_mpi = is_segmented | (target == 'deriche_mpi_baseline')
    is_omp = (target == 'deriche_omp')

    results = {}
    if is_segmented:
        # print('sw: ', sw)
        results = read_results_segments(root_dir_path)
    elif is_mpi:
        results = read_results_default(root_dir_path)
    elif is_omp:
        results = read_results_omp(root_dir_path)
    else:
        print("Invalid target parameter")
        exit(-1)

    if is_segmented:
        results = results[sw]

    # Compute arithmetic mean and variance, speedup compared to baseline, efficiency
    ws_target = []
    for i in range(len(n_cpus)):
        cpu = n_cpus[i]
        dim_results = results[dims[i]]
        dim_cpu_reps = np.array([i for i in dim_results[cpu].values()])
        dim_cpu_mean = np.mean(dim_cpu_reps)
        dim_cpu_var = np.var(dim_cpu_reps)
        ws_target.append(dim_cpu_mean)


    speedup = np.array([ b / s for (s,b) in zip(ws_target,baseline_ws )])
    print('speedup: ', speedup)

    return speedup


def efficiency_speedup_for_dim(sw, target):

    baseline_root_path = os.path.join(args.results_path, 'deriche_ref')
    baseline = read_results_ref(baseline_root_path)[dim]
    dim_baseline_reps = np.array([rep for rep in baseline.values()])
    dim_baseline_mean = np.mean(dim_baseline_reps)
    dim_baseline_var = np.var(dim_baseline_reps)

    # print('target: ', target)
    # print('baseline: ', dim_baseline_mean)

    root_dir_path = args.results_path + '/' + target
    # whether the results contains a 'SW' parameter or not
    is_segmented = (target == 'deriche_mpi_segments') | (target == 'deriche_mpi_rdma')
    is_mpi = is_segmented | (target == 'deriche_mpi_baseline')
    is_omp = (target == 'deriche_omp')

    results = {}
    if is_segmented:
        # print('sw: ', sw)
        results = read_results_segments(root_dir_path)
    elif is_mpi:
        results = read_results_default(root_dir_path)
    elif is_omp:
        results = read_results_omp(root_dir_path)
    else:
        print("Invalid target parameter")
        exit(-1)

    dim_results = {}
    if is_segmented:
        dim_results = results[sw][dim]
    else:
        dim_results = results[dim]


    # Compute arithmetic mean and variance, speedup compared to baseline, efficiency
    dim_cpus_reps = np.array([[i for i in dim_results[cpu].values()] for cpu in n_cpus])
    dim_cpus_mean = np.array([ np.mean(reps) for reps in dim_cpus_reps ])
    dim_cpus_var = np.array([ np.var(reps) for reps in dim_cpus_reps ])

    speedup = np.array([ dim_baseline_mean / runtime for runtime in dim_cpus_mean])
    efficiency = np.array([ s / cpus for (s, cpus) in zip(speedup, n_cpus)])
    # print('speedup: ', speedup)

    return speedup, efficiency

def plot_runtime(log=False):


    TITLE = 'Runtime - Strong Scaling - W,H = ' + W_POW_OF_2 +  ',' +  H_POW_OF_2

    Y_LABEL = 'Runtime [s]'

    fig, ax = plt.subplots()
    for i in range(len(targets)):
        print('plotting label', labels[i])
        runtime, variance = runtime_for_dim(sws[i], targets[i])
        plt.plot(x_ticks, runtime, line_styles[i], label=labels[i])

    ax.legend()

    # title, axes labels
    # plt.title(TITLE, y=1.07, size=12)
    plt.ylabel(' ' + Y_LABEL, rotation=0, horizontalalignment='left', y=1.02, color=AXIS_LABEL_COLOR)
    plt.xlabel(X_LABEL, color=AXIS_LABEL_COLOR)

    # remove y axis ticks
    plt.tick_params(axis='y', which='both', left=False, right=False)

    # x ticks
    ax.set_xticks(x_ticks, labels=x_labels)

    # horizontal lines, background
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.yaxis.grid(color=HORIZONTAL_LINES_COLOR, linewidth=1.0)
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    # logarithmic y axis
    if log:
        plt.semilogy(base=2)

    # hide box, add digit separator
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    filename = FILENAME
    if log:
        filename = filename + '_log'
    outfile = filename + '_ws.eps'

    print('output: ', outfile)
    plt.savefig(outfile, dpi=600, bbox_inches='tight')

def plot_speedup(log=False):

    TITLE = 'Speedup - Strong Scaling - W,H = ' + W_POW_OF_2 +  ',' +  H_POW_OF_2

    Y_LABEL = 'Speedup'

    fig, ax = plt.subplots()
    for i in range(len(targets)):
        print('plotting label=', labels[i], ' sw=', sws[i])
        if targets[i] != 'deriche_ref':
            speedup, efficiency = efficiency_speedup_for_dim(sws[i], targets[i])
            plt.plot(x_ticks, speedup, line_styles[i], label=labels[i])

    plt.plot(x_ticks, n_cpus, 'y:', label='Linear Bound')

    ax.legend(fontsize='small')

    # title, axes labels
    # plt.title(TITLE, y=1.07, size=12)
    plt.ylabel(' ' + Y_LABEL, rotation=0, horizontalalignment='left', y=1.02, color=AXIS_LABEL_COLOR)
    plt.xlabel(X_LABEL, color=AXIS_LABEL_COLOR)


    # logarithmic y axis
    if log:
        plt.semilogy(base=2)

    # remove y axis ticks
    plt.tick_params(axis='y', which='both', left=False, right=False)
    # ax.set_yticks(n_cpus)

    # x ticks
    ax.set_xticks(x_ticks, labels=x_labels)

    # horizontal lines, background
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.yaxis.grid(color=HORIZONTAL_LINES_COLOR, linewidth=1.0)
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)



    # hide box, add digit separator
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    filename = FILENAME
    if log:
        filename = filename + '_log'
    outfile = filename + '_speedup.eps'
    print('output: ', outfile)
    plt.savefig(outfile, dpi=600, bbox_inches='tight')

def plot_speedup_ws(log=False):

    TITLE = 'Speedup - Weak Scaling'

    Y_LABEL = 'Speedup'

    if dim == 14:
        dims = [1112, 12, 1213, 13, 1314, 14]
    elif dim == 13:
        dims = [11, 1112, 12, 1213, 13, 1314]
    elif dim == 12:
        dims = [1011, 11, 1112, 12, 1213, 13]
    else:
        dims = [10, 1011, 11, 1112, 12, 1213]

    fig, ax = plt.subplots()
    for i in range(len(targets)):
        print('plotting label=', labels[i], ' sw=', sws[i])
        if targets[i] != 'deriche_ref':
            speedup = speedup_for_dims(sws[i], targets[i], dims)
            plt.plot(x_ticks, speedup, line_styles[i], label=labels[i])


    plt.plot(x_ticks, n_cpus, 'y:', label='Linear Bound')

    ax.legend(fontsize='small')

    # title, axes labels
    # plt.title(TITLE, y=1.07, size=12)
    plt.ylabel(' ' + Y_LABEL, rotation=0, horizontalalignment='left', y=1.02, color=AXIS_LABEL_COLOR)
    plt.xlabel(X_LABEL, color=AXIS_LABEL_COLOR)

    # remove y axis ticks
    plt.tick_params(axis='y', which='both', left=False, right=False)

    # x ticks
    ax.set_xticks(x_ticks, labels=x_labels)

    # horizontal lines, background
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.yaxis.grid(color=HORIZONTAL_LINES_COLOR, linewidth=1.0)
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)

    # logarithmic y axis
    if log:
        plt.semilogy(base=2)

    # hide box, add digit separator
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    filename = 'deriche_' + str(dims[0])
    if log:
        filename = filename + '_log'
    outfile = filename + '_ws.eps'

    print('output: ', outfile)
    plt.savefig(outfile, dpi=600, bbox_inches='tight')


# EXECUTION
# plot_runtime()
# plot_runtime(log=True)
# plot_speedup_ws()
# plot_speedup()
plot_speedup_ws(log=True)
plot_speedup(log=True)
