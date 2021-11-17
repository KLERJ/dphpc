#!/usr/bin/env bash
# This script generates the reference we can test against
DIR=$(dirname "$0")
ABS_DIR="$(pwd)/${DIR}"
echo $ABS_DIR
ROOT_DIR="${DIR}/../.."
PB_ROOT="${ROOT_DIR}/polybench/"
datasets=('MINI_DATASET' 'SMALL_DATASET' 'MEDIUM_DATASET' 'LARGE_DATASET' 'EXTRALARGE_DATASET')
# it might be possible to set custom W and H values at compile time with -D W=1080 -D H=1920 ?
cd "$PB_ROOT"

if [ -z "$1" ]
then
  idx=0
else
  idx=$1
fi

dataset=${datasets[${idx}]}
echo $dataset
gcc -O3 -I utilities -I medley/deriche utilities/polybench.c medley/deriche/deriche.c \
 -DPOLYBENCH_DUMP_ARRAYS -D${dataset} -o deriche_ref
./deriche_ref 2> "${ABS_DIR}/deriche_$dataset.out"
