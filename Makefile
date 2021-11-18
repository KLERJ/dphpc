CFLAGS=-Wall -Wextra -Wno-comment -O3 -I ${SRC_DIR}
SIZE=MINI
BINARY_DERICHE=deriche
BINARY_SEIDEL2D_OMP=seidel-2d_omp
BINARY_HEAT3D_OMP=heat-3d_omp
SRC_DIR=./src
CC=gcc-11
MPI_CC=mpicc
EXTRA_FLAGS=-DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME -D${SIZE}_DATASET

all: deriche seidel-2d_omp heat-3d_omp

.PHONY: deriche run

deriche_mpi:
	${MPI_CC} ${CFLAGS} ${SRC_DIR}/deriche/deriche_mpi.c -DMINI_DATASET -o ${BINARY_DERICHE}

run: deriche
	mpiexec -np 2 ./${BINARY_DERICHE}

seidel-2d_omp:
	${CC} ${CFLAGS} -fopenmp -o ${BINARY_SEIDEL2D_OMP} ${SRC_DIR}/seidel-2d/seidel-2d_omp.c -I ${SRC_DIR}/seidel-2d -I polybench/utilities polybench/utilities/polybench.c ${EXTRA_FLAGS}

heat-3d_omp:
	${CC} ${CFLAGS} -fopenmp -o ${BINARY_HEAT3D_OMP} ${SRC_DIR}/heat-3d/heat-3d_omp.c -I ${SRC_DIR}/heat-3d -I polybench/utilities polybench/utilities/polybench.c ${EXTRA_FLAGS}

clean:
	@ rm -f ${BINARY_DERICHE} ${BINARY_SEIDEL2D_OMP} ${BINARY_HEAT3D_OMP}
