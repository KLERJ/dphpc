#!/usr/bin/env bash
set -e

if [ "$#" -ne 2 ]; then
  echo "Usage: ./benchmark_deriche.sh <n-reps> <result-dir>"
  exit 1
fi
N_REPS=$1
RESULT_DIR=$2

# DATASETS=('MINI' 'SMALL' 'MEDIUM' 'LARGE' 'EXTRALARGE')
DERICHE_DIMS=(100 1000 10000 20000 40000)

#######################################
# Perform benchmark on the target.
# Globals:
#   DERICHE_DIMS
#   N_REPS
#   RESULT_DIR
# Arguments:
#   The target to benchmark
# Outputs:
#   Writes result in $target.$size.$rep
#######################################
benchmark() {
  target=$1
  # for size in "${DATASETS[@]}"; do
  for dim in ${DERICHE_DIMS[@]}; do
    make clean
    make $target DERICHE_DIM=$dim
    for rep in $(seq 1 $N_REPS); do
      bin/$target > $RESULT_DIR/$target.$dim.$rep
    done
  done
}


set -x
mkdir -p $RESULT_DIR

# benchmark deriche_ref
benchmark deriche_omp
