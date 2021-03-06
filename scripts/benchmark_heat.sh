#!/bin/bash

#Problem setup

N_THREADS=(1 2 4 8 16 32 64 128 216 256)

MPI_OMP_THREADS=(4 8 16 32 64 128 256)
MPI_THREADS=(2 2 2 4 4 8 16)
OMP_THREADS=(2 4 8 8 16 16 16)

RUN_DIR=$1     # Root directory of the run
USE_SCOREP=$2  # Use Score-P true/false

DN=$3  		   # Problem Size
DSTEPS=$4      # Number of Steps




function run_heat3d_mpi(){


	target=$1
	use_scorep=$2
	extra_flags=$3
	skip_first=${4:-false}

	if [ "$skip_first" = "true" ]; then
		start_zero=1
	else
		start_zero=0
	fi

	if [ "$use_scorep" = "true" ]; then
		echo "Using Scorep"
		CC='scorep mpicc'
	else
		CC=mpicc
	fi
	
	# make clean
	make $target MPI_CC="$CC" EXTRA_FLAGS="-DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME $extra_flags" BIN_DIR=$RUN_DIR;
	
	EXEC=$(printf "%s_%s" $target ${extra_flags: -1})

	mv $RUN_DIR/$target $RUN_DIR/$EXEC
	
	hwloc-ls --output-format xml > $RUN_DIR/hwloc.xml 
	lscpu > $RUN_DIR/cpu.txt

	for ((i=$start_zero;i<${#N_THREADS[@]};i++))
	do
		n_threads=${N_THREADS[$i]}
	
		EXP_NAME=$(printf "%s_%dx_%dsize%s" $target $n_threads $DN $extra_flags)

	
		EXP_DIR="$RUN_DIR/$EXP_NAME"
		mkdir -p $EXP_DIR
		mkdir $EXP_DIR/trace
	

		export SCOREP_ENABLE_PROFILING=true
		export SCOREP_ENABLE_TRACING=true
		export SCOREP_EXPERIMENT_DIRECTORY="$EXP_DIR/trace"
	
		echo "Job $i $target running on $n_threads cores"

		time mpirun -np $n_threads $RUN_DIR/$EXEC $DN $DSTEPS $EXP_DIR/output_$DN.bin $EXP_DIR > $EXP_DIR/output.txt 2> $EXP_DIR/error.txt

	done
}





function run_heat3d_mpi_omp_hybrid(){


	target=$1
	use_scorep=$2
	extra_flags=$3
	

	if [ "$use_scorep" = "true" ]; then
		echo "Using Scorep"
		CC='scorep mpicc'
	else
		CC=mpicc
	fi
	
	# make clean
	make $target MPI_CC="$CC" EXTRA_FLAGS=" -DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME $extra_flags" BIN_DIR=$RUN_DIR;
	hwloc-ls --output-format xml > $RUN_DIR/hwloc.xml 
	lscpu > $RUN_DIR/cpu.txt

	for ((i=0;i<${#MPI_OMP_THREADS[@]};i++))
	do
		n_threads=${MPI_OMP_THREADS[$i]}
		mpi_threads=${MPI_THREADS[$i]}
		omp_threads=${OMP_THREADS[$i]}
	
		EXP_NAME=$(printf "%s_%dx_%dsize_%dx%d_%s" $target $n_threads $DN $mpi_threads $omp_threads $extra_flags)

	
		EXP_DIR="$RUN_DIR/$EXP_NAME"
		mkdir -p $EXP_DIR
		mkdir $EXP_DIR/trace
	

		export SCOREP_ENABLE_PROFILING=true
		export SCOREP_ENABLE_TRACING=true
		export SCOREP_EXPERIMENT_DIRECTORY="$EXP_DIR/trace"
		export OMP_NUM_THREADS=$omp_threads
	
		echo "Job $i $target running on $mpi_threads x $omp_threads = $n_threads cores"
		
		unset LSB_AFFINITY_HOSTFILE
		# time mpirun -n $mpi_threads --map-by node:PE=$omp_threads $RUN_DIR/$target $DN $DSTEPS $EXP_DIR/output_$DN.bin $EXP_DIR > $EXP_DIR/output.txt 2> $EXP_DIR/error.txt
		time mpirun -n $mpi_threads --loadbalance --cpus-per-proc $omp_threads $RUN_DIR/$target $DN $DSTEPS $EXP_DIR/output_$DN.bin $EXP_DIR > $EXP_DIR/output.txt 2> $EXP_DIR/error.txt

	done
}


function run_heat3d_omp(){


	target=$1

	skip_first=${2:-false}

	if [ "$skip_first" = "true" ]; then
		start_zero=1
	else
		start_zero=0
	fi

	
	# make clean
	hwloc-ls --output-format xml > $RUN_DIR/hwloc.xml 
	lscpu > $RUN_DIR/cpu.txt

	for ((i=0;i<${#N_THREADS[@]};i++))
	do
		n_threads=${N_THREADS[$i]}

		export OMP_NUM_THREADS=$n_threads
		export MKL_NUM_THREADS=$n_threads
	
		EXP_NAME=$(printf "%s_%dx_%dsize" $target $n_threads $DN)

	
		EXP_DIR="$RUN_DIR/$EXP_NAME"
		mkdir -p $EXP_DIR
		mkdir $EXP_DIR/trace
	
	    make $target EXTRA_FLAGS="-DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME -DN=$DN -DTSTEPS=$DSTEPS" BIN_DIR=$EXP_DIR; 

	
		echo "Job $i $target running on $n_threads cores"

		$EXP_DIR/$target $DN $DSTEPS $EXP_DIR/output_$DN.bin $EXP_DIR > $EXP_DIR/output.txt 2> $EXP_DIR/error.txt

	done
}


function run_heat3d_baseline(){
	cd polybench/stencils/heat-3d

	make clean
	make EXTRA_FLAGS=" -DPOLYBENCH_TIME -DN=$DN -DTSTEPS=$DSTEPS" BIN_DIR=$RUN_DIR

	time $RUN_DIR/heat-3d > $RUN_DIR/baseline_output.txt 2> $RUN_DIR/baseline_error.txt
	cd ../../..
}

run_heat3d_mpi_omp_hybrid heat-3d_mpi_omp $USE_SCOREP ""

run_heat3d_mpi heat-3d_mpi $USE_SCOREP "-DZ_DIM=0" false 
run_heat3d_mpi heat-3d_mpi $USE_SCOREP "-DZ_DIM=1" true # Skip 1 Core iter
run_heat3d_mpi heat-3d_mpi $USE_SCOREP "-DZ_DIM=2" true # Skip 1 Core iter


# run_heat3d_mpi heat-3d_mpi_avx2 $USE_SCOREP "" true # Skip 1 Core iter
run_heat3d_baseline;

run_heat3d_omp heat-3d_omp true; # Skip 1 core iter
