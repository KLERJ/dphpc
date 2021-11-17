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
gcc -O3 -I utilities -I stencils/seidel-2d utilities/polybench.c stencils/seidel-2d/seidel-2d.c \
 -DPOLYBENCH_DUMP_ARRAYS -D${dataset} -o seidel-2d_ref
./seidel-2d_ref 2> ${ABS_DIR}/seidel-2d_$dataset.out
