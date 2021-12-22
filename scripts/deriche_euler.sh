#!/usr/bin/env bash
set -e
set -x
SCRIPTS_DIR=$(dirname "$0")
cd $SCRIPTS_DIR
pwd

bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 24:00 -J "deriche" -n 32 <<EOF
./benchmark_deriche.sh 50 deriche_mpi_baseline
./benchmark_deriche.sh 50 deriche_mpi_rdma
./benchmark_deriche.sh 50 deriche_mpi_segments
./benchmark_deriche.sh 50 deriche_mpi_avx
./benchmark_deriche.sh 50 deriche_omp
./benchmark_deriche.sh 50 deriche_ref
EOF
