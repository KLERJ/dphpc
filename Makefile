CFLAGS=-Wall -Wextra -Wno-comment -O3 -I ${SRC_DIR} -D${SIZE}_DATASET
SIZE=MINI
BINARY_DERICHE=deriche
SRC_DIR=./src
CC=gcc
MPI_CC=mpicc

.PHONY: deriche run

deriche:
	${MPI_CC} ${CFLAGS} ${SRC_DIR}/deriche.c -DMINI_DATASET -o ${BINARY_DERICHE}

run: deriche
	mpiexec -np 2 ./${BINARY_DERICHE}
