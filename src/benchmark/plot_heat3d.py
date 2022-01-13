import os
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

# Set these for your machine:
rootdir = '/Users/lukevolpatti/Downloads/s'
outdir  = '/Users/lukevolpatti/Downloads/plots'

# Parameters
programs = ['mpi', 'omp']

# Graph properties
AXIS_LABEL_COLOR = (.4, .4, .4)
BACKGROUND_COLOR = (.88, .88, .88)
HORIZONTAL_LINES_COLOR = (1, 1, 1)


def plot_by_size(subfolder):
    subfolder_info = subfolder.split("_")
    num_iters = subfolder_info[-2][:-1]
    prob_size = subfolder_info[-3]

    mpid0_total_runtime = {}
    mpid1_total_runtime = {}
    mpid2_total_runtime = {}
    omp_total_runtime = {}
    mpid0_per_iter = {}
    mpid1_per_iter = {}
    mpid2_per_iter = {}
    omp_per_iter = {}
    mpid0_comm_time = {}
    mpid1_comm_time = {}
    mpid2_comm_time = {}
    baseline_runtime = 0
    dim = -1


    instance_dirs = [ f.path for f in os.scandir(subfolder) if f.is_dir() ]

    for instance_dir in instance_dirs:
        instance_dir_info = instance_dir.split('_')
        if 'DIM' in instance_dir_info[-1] or 'x' in instance_dir_info[-1]:
            nprocs = int(instance_dir_info[-3][:-1])
            program_type = instance_dir_info[-4]
            if 'DIM=0' in instance_dir_info[-1]:
                dim = 0
            elif 'DIM=1' in instance_dir_info[-1]:
                dim = 1
            elif 'DIM=2' in instance_dir_info[-1]:
                dim = 2
        else:
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
                if dim == 0:
                    mpid0_total_runtime[nprocs] = float(total_runtime)
                elif dim == 1:
                    mpid1_total_runtime[nprocs] = float(total_runtime)
                elif dim == 2:
                    mpid2_total_runtime[nprocs] = float(total_runtime)
            elif program_type == "omp":
                omp_total_runtime[nprocs] = float(total_runtime)
            elif program_type == 'baseline':
                baseline_runtime = float(total_runtime)


        # Get per processor per iteration runtimes
        benchmark_files = [ f.path for f in os.scandir(instance_dir) if f.is_file() ]
        for benchmark_file in benchmark_files:
            if 'benchmark_' in benchmark_file:
                with open(benchmark_file, 'r') as f:
                    if program_type == 'mpi':
                        lines = f.read().splitlines()
                        values = lines[-1].split(',')
                        values = np.array([float(i) for i in values])
                        comm_comp_time = lines[2].split(',')
                        comm_comp_time = np.array([float(i) for i in comm_comp_time])
                        comp_time = lines[3].split(',')
                        comp_time = np.array([float(i) for i in comp_time])
                        comm_time = comm_comp_time - comp_time
                        
                        if dim == 0:
                            if nprocs in mpid0_per_iter:
                                mpid0_per_iter[nprocs] = np.append(mpid0_per_iter[nprocs], values)
                            else:
                                mpid0_per_iter[nprocs] = values
                            if nprocs in mpid0_comm_time:
                                mpid0_comm_time[nprocs] = np.append(mpid0_comm_time[nprocs], comm_time)
                            else:
                                mpid0_comm_time[nprocs] = values
                        if dim == 1:
                            if nprocs in mpid1_per_iter:
                                mpid1_per_iter[nprocs] = np.append(mpid1_per_iter[nprocs], values)
                            else:
                                mpid1_per_iter[nprocs] = values
                            if nprocs in mpid1_comm_time:
                                mpid1_comm_time[nprocs] = np.append(mpid1_comm_time[nprocs], comm_time)
                            else:
                                mpid1_comm_time[nprocs] = values
                        if dim == 2:
                            if nprocs in mpid2_per_iter:
                                mpid2_per_iter[nprocs] = np.append(mpid2_per_iter[nprocs], values)
                            else:
                                mpid2_per_iter[nprocs] = values
                            if nprocs in mpid2_comm_time:
                                mpid2_comm_time[nprocs] = np.append(mpid2_comm_time[nprocs], comm_time)
                            else:
                                mpid2_comm_time[nprocs] = values
                    elif program_type == 'omp':
                        values = f.read().splitlines()[0].split(',')
                        values = np.array([float(i) for i in values])
                        values = values[2:]
                        if nprocs in omp_per_iter:
                            omp_per_iter[nprocs] = np.append(omp_per_iter[nprocs], values)
                        else:
                            omp_per_iter[nprocs] = values
        
    print(mpid0_total_runtime)
    print(mpid1_total_runtime)
    print(mpid2_total_runtime)
    print(omp_total_runtime)

    # Total runtime plots
    mpi_l = sorted(mpid0_total_runtime.items())
    mpi_x, mpi_y = zip(*mpi_l)
    if omp_total_runtime:
        omp_l = sorted(omp_total_runtime.items())
        omp_x, omp_y = zip(*omp_l)

    fig, ax = plt.subplots()
    plt.plot(mpi_x, mpi_y, 'b.-', label='mpi')
    if omp_total_runtime:
        plt.plot(omp_x, omp_y, 'g.-', label='omp')
    plt.legend()
    plt.xlabel("# cores", color=AXIS_LABEL_COLOR)
    plt.ylabel("Total runtime (s)", color=AXIS_LABEL_COLOR)
    plt.title("heat-3d, problem size: " + prob_size + ", iterations: " + num_iters)

    plt.tick_params(axis='y', which='both', left=False, right=False)
    plt.grid(axis='y', color=HORIZONTAL_LINES_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    plt.savefig("heat-3d_total_runtime_" + prob_size + "_" + num_iters + "x.png")

    # Communication time histograms
    for procs in mpid0_comm_time:
        if procs != 1:
            nbins = 15
            times0 = mpid0_comm_time[procs]
            times1 = mpid1_comm_time[procs]
            times2 = mpid2_comm_time[procs]
            fig, ax = plt.subplots()
            plt.hist([times0, times1, times2], bins=nbins, histtype='bar', label=['d0','d1','d2'])
            plt.legend()
            plt.xlabel('Time (s)')
            plt.title('Problem size: ' + prob_size + ', cores: ' + str(procs) + ', iterations: ' + num_iters)
            #plt.hist(times1, bins=nbins, histtype='bar', stacked=True)
            #plt.hist(times2, bins=nbins, histtype='bar', stacked=True)
            plt.savefig("heat-3d_comm_time_nprocs" + str(procs) + "_" + prob_size + "_" + num_iters + "x.png")

    # Speedup plot
    if '1024' in subfolder:
        mpid0_l = sorted(mpid0_total_runtime.items())
        mpid0_x, mpid0_y = zip(*mpid0_l)
        mpid0_y = np.array(mpid0_y)
        mpid0_y = baseline_runtime / mpid0_y

        mpid1_l = sorted(mpid1_total_runtime.items())
        mpid1_x, mpid1_y = zip(*mpid1_l)
        mpid1_y = np.array(mpid1_y)
        mpid1_y = baseline_runtime / mpid1_y

        mpid2_l = sorted(mpid2_total_runtime.items())
        mpid2_x, mpid2_y = zip(*mpid2_l)
        mpid2_y = np.array(mpid2_y)
        mpid2_y = baseline_runtime / mpid2_y

        if omp_total_runtime:
            omp_l = sorted(omp_total_runtime.items())
            omp_x, omp_y = zip(*omp_l)
            omp_y = np.array(omp_y)
            omp_y = baseline_runtime / omp_y

        fig, ax = plt.subplots()
        plt.plot(mpid0_x, mpid0_y, 'b.-', label='mpid0')
        plt.plot(mpid1_x, mpid1_y, 'y.-', label='mpid1')
        plt.plot(mpid2_x, mpid2_y, 'r.-', label='mpid2')
        if omp_total_runtime:
            plt.plot(omp_x, omp_y, 'g.-', label='omp')
        plt.legend()
        plt.xlabel("# cores", color=AXIS_LABEL_COLOR)
        plt.ylabel("Speedup", color=AXIS_LABEL_COLOR)
        plt.title("heat-3d, problem size: " + prob_size + ", iterations: " + num_iters)

        plt.tick_params(axis='y', which='both', left=False, right=False)
        plt.grid(axis='y', color=HORIZONTAL_LINES_COLOR)
        ax.set_facecolor(BACKGROUND_COLOR)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

        plt.savefig("heat-3d_speedup_" + prob_size + "_" + num_iters + "x.png")

    
    # Per processor per iteration plots
    mpi_l = sorted(mpid0_per_iter.items())
    mpi_x, mpi_y = zip(*mpi_l)
    mpi_err = [np.std(i) for i in mpi_y]
    mpi_y = [np.mean(i) for i in mpi_y]

    if omp_per_iter:
        omp_l = sorted(omp_per_iter.items())
        omp_x, omp_y = zip(*omp_l)
        omp_err = [np.std(i) for i in omp_y]
        omp_y = [np.mean(i) for i in omp_y]

    fig, ax = plt.subplots()
    plt.errorbar(mpi_x, mpi_y,
            yerr=mpi_err, 
            ecolor='black', 
            fmt = 'b.-',
            label='mpi', 
            capsize=2)
    if omp_per_iter:
        plt.errorbar(omp_x, omp_y,
                yerr=omp_err, 
                ecolor='black', 
                fmt = 'g.-',
                label='omp', 
                capsize=2)
    #plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel("# cores", color=AXIS_LABEL_COLOR)
    plt.ylabel("Time per iteration per core (s)", color=AXIS_LABEL_COLOR)
    plt.title("heat-3d, problem size: " + prob_size + ", iterations: " + num_iters)
    plt.tick_params(axis='y', which='both', left=False, right=False)
    plt.grid(axis='y', color=HORIZONTAL_LINES_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    plt.savefig("heat-3d_per_iter_" + prob_size + "_" + num_iters + "x.png")

    return mpid0_total_runtime

    
if __name__ == '__main__':
    subfolders = [ f.path for f in os.scandir(rootdir) if f.is_dir() ]
    print(subfolders)

    total_runtimes = {}
    problem_sizes = [1024, 1536, 2048]
    nprocs = [8, 16, 32,64, 128, 216, 256]

    for subfolder in subfolders:
        times = plot_by_size(subfolder)
        if '1024' in subfolder:
            total_runtimes[1024] = times
        if '1536' in subfolder:
            total_runtimes[1536] = times
        if '2048' in subfolder:
            total_runtimes[2048] = times
    print(total_runtimes)

    # Making weak scaling plots
    fig, ax = plt.subplots()
    for procs in nprocs:
        time_by_proc = []
        for psize in problem_sizes:
            time_by_proc.append(total_runtimes[psize][procs])
        plt.plot(problem_sizes, time_by_proc, '.-', label=str(procs)+' procs')

    plt.legend()
    plt.xlabel('Problem size')
    plt.ylabel('Total runtime (s)')
    plt.savefig("weak_scaling.png")
