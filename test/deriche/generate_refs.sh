#!/usr/bin/env bash
# This scripts generates the reference we can test against
DIR=$(dirname "$0")
ABS_DIR="$(pwd)/${DIR}"
echo $ABS_DIR
ROOT_DIR="${DIR}/../.."
PB_ROOT="${ROOT_DIR}/polybench/"
datasets=('MINI_DATASET' 'SMALL_DATASET' 'MEDIUM_DATASET' 'LARGE_DATASET' 'EXTRALARGE_DATASET')
# it might be possible to set custome W and H values at compile time with -D W=1080 -D H=1920 ?
cd $PB_ROOT


dataset=${datasets[0]}
echo $dataset
gcc -O3 -I utilities -I medley/deriche utilities/polybench.c medley/deriche/deriche.c \
 -DPOLYBENCH_DUMP_ARRAYS -D${dataset} -o deriche_ref
./deriche_ref 2> ${ABS_DIR}/deriche_$dataset.out
