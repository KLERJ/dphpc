SIZE=SMALL
CFLAGS=-Wall -Wextra -Wno-comment -O3 -I polybench/utilities polybench/utilities/polybench.c
DUMP=
EXTRA_FLAGS=-DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME -D${SIZE}_DATASET ${DUMP}

SRC_DIR=./src
BIN_DIR=./bin
BINARY_DERICHE=${BIN_DIR}/deriche
BINARY_DERICHE_OMP=${BIN_DIR}/deriche_omp
BINARY_SEIDEL2D_OMP=${BIN_DIR}/seidel-2d_omp
BINARY_HEAT3D_OMP=${BIN_DIR}/heat-3d_omp
BINARY_HEAT3D_MPI=${BIN_DIR}/heat-3d_mpi

$(shell mkdir -p bin)
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CC=gcc
endif
ifeq ($(UNAME), Darwin)
CC=gcc-11
endif

MPI_CC=mpicc

all: deriche deriche_omp seidel-2d_omp heat-3d_omp heat-3d_mpi

.PHONY: clean run

DERICHE_DIM=1000
deriche:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS}  -o ${BINARY_DERICHE} ${SRC_DIR}/deriche/deriche_mpi.c

deriche_omp:
	${CC} ${CFLAGS} ${EXTRA_FLAGS}  -fopenmp -o ${BINARY_DERICHE_OMP} ${SRC_DIR}/deriche/deriche_omp.c -DW=${DERICHE_DIM} -DH=${DERICHE_DIM}

deriche_ref:
	${CC} ${CFLAGS} ${EXTRA_FLAGS}  -o ${BIN_DIR}/deriche_ref polybench/medley/deriche/deriche.c -DW=${DERICHE_DIM} -DH=${DERICHE_DIM}

run: deriche
	mpiexec -np 2 ./${BINARY_DERICHE}

seidel-2d_omp:
	${CC} ${CFLAGS}  -fopenmp -o ${BINARY_SEIDEL2D_OMP} ${SRC_DIR}/seidel-2d/seidel-2d_omp.c -I ${SRC_DIR}/seidel-2d ${EXTRA_FLAGS}

heat-3d_omp:
	${CC} ${CFLAGS}  -fopenmp -o ${BINARY_HEAT3D_OMP} ${SRC_DIR}/heat-3d/heat-3d_omp.c -I ${SRC_DIR}/heat-3d ${EXTRA_FLAGS}

heat-3d_mpi:
	${MPI_CC} ${CFLAGS}  ${SRC_DIR}/heat-3d/heat-3d_mpi.c  -o ${BINARY_HEAT3D_MPI} ${EXTRA_FLAGS}

clean:
	@ rm -f ${BINARY_DERICHE} ${BINARY_DERICHE_OMP} ${BINARY_SEIDEL2D_OMP} ${BINARY_HEAT3D_OMP}
