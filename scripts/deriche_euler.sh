#!/usr/bin/env bash
set -e
set -x
SCRIPTS_DIR=$(dirname "$0")
cd $SCRIPTS_DIR
pwd

# Baseline deriche.
bsub -R "model=XeonE3_1585Lv5 rusage[mem=4000]" -W 4:00 -J "deriche-ref" -n 32 ./benchmark_deriche.sh 5 deriche_mpi_rdma
