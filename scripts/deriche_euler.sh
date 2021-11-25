#!/usr/bin/env bash
set -e
set -x
SCRIPTS_DIR=$(dirname "$0")
cd $SCRIPTS_DIR
pwd

# Baseline deriche.
bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J "deriche-ref" -n 1 ./benchmark_deriche.sh 5 deriche_ref deriche-ref

for ncpus in {1,2,4,8,16,32}; do
  job_name="deriche-omp-$ncpus"
  bsub -R "model=EPYC_7H12 rusage[mem=4000]" -W 4:00 -J $job_name -n $ncpus ./benchmark_deriche.sh 5 deriche_omp $job_name
done
