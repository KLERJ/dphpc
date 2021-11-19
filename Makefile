SIZE=LARGE
CFLAGS=-Wall -Wextra -Wno-comment -O3 -I polybench/utilities
EXTRA_FLAGS=-DPOLYBENCH_USE_C99_PROTO -DPOLYBENCH_TIME -D${SIZE}_DATASET

SRC_DIR=./src
BIN_DIR=./bin
BINARY_DERICHE=${BIN_DIR}/deriche
BINARY_DERICHE_OMP=${BIN_DIR}/deriche_omp
BINARY_SEIDEL2D_OMP=${BIN_DIR}/seidel-2d_omp
BINARY_HEAT3D_OMP=${BIN_DIR}/heat-3d_omp

$(shell mkdir -p bin)
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CC=gcc
endif
ifeq ($(UNAME), Darwin)
CC=gcc-11
endif

MPI_CC=mpicc

OBJS=polybench.o

all: deriche deriche_omp seidel-2d_omp heat-3d_omp

.PHONY: clean run

deriche: ${OBJS}
	${MPI_CC} ${CFLAGS} ${OBJS} ${SRC_DIR}/deriche/deriche_mpi.c -DMINI_DATASET -o ${BINARY_DERICHE}

deriche_omp: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -fopenmp -o ${BINARY_DERICHE_OMP}  ${SRC_DIR}/deriche/deriche_omp.c

run: deriche
	mpiexec -np 2 ./${BINARY_DERICHE}

seidel-2d_omp: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -fopenmp -o ${BINARY_SEIDEL2D_OMP} ${SRC_DIR}/seidel-2d/seidel-2d_omp.c -I ${SRC_DIR}/seidel-2d ${EXTRA_FLAGS}

heat-3d_omp: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -fopenmp -o ${BINARY_HEAT3D_OMP} ${SRC_DIR}/heat-3d/heat-3d_omp.c -I ${SRC_DIR}/heat-3d ${EXTRA_FLAGS}

polybench.o:
	${CC} ${CFLAGS} ${EXTRA_FLAGS} -c polybench/utilities/polybench.c

clean:
	@ rm -f ${OBJS} ${BINARY_DERICHE} ${BINARY_DERICHE_OMP} ${BINARY_SEIDEL2D_OMP} ${BINARY_HEAT3D_OMP}
