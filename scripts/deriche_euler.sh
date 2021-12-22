#!/usr/bin/env bash
set -e
set -x
SCRIPTS_DIR=$(dirname "$0")
cd $SCRIPTS_DIR
pwd

# bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_mpi_baseline" -n 32 ./benchmark_deriche.sh 25 deriche_mpi_baseline
# bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_mpi_rdma" -n 32 ./benchmark_deriche.sh 25 deriche_mpi_rdma
bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_mpi_segments" -n 32 ./benchmark_deriche.sh 25 deriche_mpi_segments
bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_mpi_avx" -n 32 ./benchmark_deriche.sh 25 deriche_mpi_avx
bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_omp" -n 32 ./benchmark_deriche.sh 25 deriche_omp
bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche_ref" -n 1 ./benchmark_deriche.sh 25 deriche_ref
