#!/usr/bin/env bash
set -e

if [ "$#" -ne 3 ]; then
  echo "Usage: ./benchmark_deriche.sh <n-reps> <target> <job-name>"
  exit 1
fi
N_REPS=$1
TARGET=$2
JOB_NAME=$3

# DATASETS=('MINI' 'SMALL' 'MEDIUM' 'LARGE' 'EXTRALARGE')
DERICHE_DIMS=(100 1000 10000 20000 40000)

set -x
mkdir -p $JOB_NAME
# for size in "${DATASETS[@]}"; do
for dim in ${DERICHE_DIMS[@]}; do
  make clean
  make $TARGET DERICHE_DIM=$dim BINARY_DERICHE_OMP=./bin/$JOB_NAME
  for rep in $(seq 1 $N_REPS); do
    bin/$TARGET > $JOB_NAME/$TARGET.$dim.$rep
  done
done



