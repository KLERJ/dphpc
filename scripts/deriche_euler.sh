#!/usr/bin/env bash
set -e
set -x
SCRIPTS_DIR=$(dirname "$0")
cd $SCRIPTS_DIR
pwd

# Baseline deriche.
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 24:00 -J "deriche_mpi_baseline" -n 32 ./benchmark_deriche.sh 50 deriche_mpi_baseline
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 24:00 -J "deriche_mpi_rdma" -n 32 ./benchmark_deriche.sh 50 deriche_mpi_rdma
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 24:00 -J "deriche_mpi_segments" -n 32 ./benchmark_deriche.sh 50 deriche_mpi_segments
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 24:00 -J "deriche_omp" -n 32 ./benchmark_deriche.sh 50 deriche_omp
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 24:00 -J "deriche_ref" -n 1 ./benchmark_deriche.sh 50 deriche_ref
