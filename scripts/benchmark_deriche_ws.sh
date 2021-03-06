#!/usr/bin/env bash
# set -e

# Weak scaling experiments for deriche - H = 2W

if (($# != 2 )); then
  echo "Usage: ./benchmark_deriche_ws.sh <n-reps> <target>"
  exit 1
fi
N_REPS=$1
TARGET=$2

# Assume using the lmod system
module load gcc/8.2.0 openmpi/4.0.2

# DATASETS=('MINI' 'SMALL' 'MEDIUM' 'LARGE' 'EXTRALARGE')
DERICHE_DIMS=(10 11 12 13)

SCRIPTS_DIR="$(dirname "$0")"
PROJECT_DIR="${SCRIPTS_DIR}/.."

set -x
cd $PROJECT_DIR
pwd

results_base="results/$TARGET"

function benchmark_ref(){
  for dim in ${DERICHE_DIMS[@]}; do
    make clean
    dim2=$((dim + 1))
    make $TARGET DIM_W=$((2 ** $dim)) DIM_H=$((2 ** $dim2))

    for rep in $(seq 1 $N_REPS); do
        results_dir=$results_base/dim_$dim$dim2/rep_$rep
        mkdir -p $results_dir
        bin/$TARGET > $results_dir/stdout 2>&1
    done
  done
}

function benchmark() {
  for dim in "${DERICHE_DIMS[@]}"; do
    make clean
    dim2=$((dim + 1))
    make $TARGET DIM_W=$((2 ** $dim)) DIM_H=$((2 ** $dim2))

    for ncpus in {1,2,4,8,16,32}; do

      for rep in $(seq 1 "${N_REPS}"); do
        results_dir="${results_base}/dim_${dim}${dim2}/cpu_${ncpus}/rep_${rep}"
        mkdir -p "${results_dir}"
        if [[ "${TARGET}" =~ .*"omp".* ]]; then
          export OMP_NUM_THREADS="${ncpus}"
          bin/"${TARGET}" > "${results_dir}/stdout"
        elif [[ "${TARGET}" =~ .*"mpi".* ]]; then
          mpirun -n "${ncpus}" "bin/${TARGET}" "${results_dir}" > "$results_dir/stdout" 2>&1
        fi
      done
    done
  done
}

function benchmark_segments() {

  SWs=(1 2 4 8 16 32)
  for sw in ${SWs[@]}; do
    for dim in ${DERICHE_DIMS[@]}; do
        make clean
        dim2=$((dim + 1))
        make $TARGET DIM_W=$((2 ** $dim)) DIM_H=$((2 ** $dim2)) SEG_WIDTH=${sw}

          for ncpus in {1,2,4,8,16,32}; do
            for rep in $(seq 1 $N_REPS); do
              results_dir=$results_base/sw_$sw/dim_$dim$dim2/cpu_$ncpus/rep_$rep
              mkdir -p $results_dir
              mpirun -n $ncpus bin/$TARGET $results_dir > $results_dir/stdout 2>&1
            done
          done
      done
  done
}

if [[ $TARGET == deriche_ref ]]; then
  benchmark_ref
elif [[ $TARGET == deriche_mpi_segments || $TARGET == deriche_mpi_rdma ]]; then
  benchmark_segments
else
  benchmark
fi
