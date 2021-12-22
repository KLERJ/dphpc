#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse
import os

########### GRAPH PROPERTIES ############

TITLE = 'Runtime deriche. LARGE dataset. Apple M1 Pro, avg. of 10 runs'
Y_LABEL = 'Runtime [seconds]'
X_LABEL = '# cores'
FILENAME = 'test'
AXIS_LABEL_COLOR = (.4, .4, .4)
BACKGROUND_COLOR = (.88, .88, .88)
HORIZONTAL_LINES_COLOR = (1, 1, 1)

# ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument('results_path', help='root directory of the benchmark results')
parser.add_argument('target', help='target version')
args = parser.parse_args()

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
                    bm_results.append(bm_file_res)
                # Add results of this rep to results of this cpu number
                cpu_results[rep] = bm_results
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
                        file = open(bm_file_path, 'r')
                        for line in file:
                            bm_file_res.append(float(line))
                        bm_results.append(bm_file_res)
                    # Add results of this rep to results of this cpu number
                    cpu_results[rep] = bm_results
                # Add results of this number of cpu to results of this dim
                dim_results[cpus] = cpu_results
            sw_results[dim] = dim_results
        results[sw] = sw_results

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
            print(bm_file_path)
            print(bm_results)
            dim_results[rep] = bm_results
        results[dim] = dim_results
    return results

FILENAME = args.target + "_plot"
# Extract data from files in results dir
root_dir_path = args.results_path + '/' + args.target
results = {}
if (args.target == 'deriche_mpi_segments') | (args.target == 'deriche_mpi_'):
    results = read_results_segments(root_dir_path)
elif args.target == 'deriche_mpi_baseline':
    results = read_results_default(root_dir_path)
elif args.target == 'deriche_ref':
    results = read_results_ref(root_dir_path)
else:
    print("Invalid target parameter")
    exit(-1)

print(results)



procs = ['1', '2', '4', '8', '16', '32']
reference_impl = [0.1356376 for i in range(len(procs))]
baseline = [0.21723179999999997, 0.09163249999999999, 0.0488108, 0.033523599999999994]
rdma = [0.2195272, 0.0859053, 0.044615699999999994, 0.0324671]
data = [reference_impl, baseline, rdma]

#########################################

procs = [int(x) for x in procs]
line_styles = ['bo-', 'ro-', 'go-', 'yo-', 'mo-']

fig, ax = plt.subplots()
for i in range(len(data)):
    plt.plot(procs, data[i], line_styles[i])

# line labels
plt.text(x=5, y=0.14, s='ref impl', color='blue')
plt.text(x=5, y=0.05, s='MPI baseline', color='red')
plt.text(x=3, y=0.03, s='MPI RDMA', color='green')

# axes
# plt.ylim(0.025, 0.25)
# plt.xlim(0, 9)
# plt.yticks(np.arange(0, 2.76, .25))

# title, axes labels
plt.title(TITLE, y=1.07, loc='left', x=-0.1)
plt.ylabel(' ' + Y_LABEL, rotation=0, horizontalalignment='left', y=1.02, color=AXIS_LABEL_COLOR)
plt.xlabel(X_LABEL, color=AXIS_LABEL_COLOR)

# remove y axis ticks
plt.tick_params(axis='y', which='both', left=False, right=False)

# horizontal lines, background
plt.grid(axis='y', color=HORIZONTAL_LINES_COLOR)
ax.set_facecolor(BACKGROUND_COLOR)

# hide box, add digit separator
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

# save
# plt.savefig(FILENAME + '.svg')
plt.savefig(FILENAME + '.png', dpi=300)
