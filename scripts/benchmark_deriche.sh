#!/usr/bin/env bash
set -e

if (($# != 2 )); then
  echo "Usage: ./benchmark_deriche.sh <n-reps> <target>"
  exit 1
fi
N_REPS=$1
TARGET=$2

# DATASETS=('MINI' 'SMALL' 'MEDIUM' 'LARGE' 'EXTRALARGE')
DERICHE_DIMS=(100 1000 10000 20000 40000)

SCRIPTS_DIR="$(dirname "$0")"
PROJECT_DIR="${SCRIPTS_DIR}/.."

set -x
cd $PROJECT_DIR
pwd

results_base="results/$TARGET"

# for size in "${DATASETS[@]}"; do
for dim in "${DERICHE_DIMS[@]}"; do
  make clean
  make "${TARGET}" DERICHE_DIM="${dim}"

  for ncpus in {1,2,4,8,16,32}; do

    for rep in "$(seq 1 "${N_REPS}")"; do
      results_dir="${results_base}/dim-${dim}/cpu-${ncpus}/rep-${rep}"
      mkdir -p "${results_dir}"

      if [[ "${TARGET}" =~ .*"omp".* ]]; then
        export OMP_NUM_THREADS="${ncpus}"
        bin/"${TARGET}" > "${results_dir}/stdout"
      elif [[ "${TARGET}" =~ .*"mpi".* ]]; then
        mpirun -n "${ncpus}" "bin/${TARGET}" "${results_dir}" > "$results_dir/stdout"
      fi
    done

  done
done



