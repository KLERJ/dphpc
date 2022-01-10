SIZE=MINI
CFLAGS=-Wall -Wextra -Wno-comment -O3 -march=native -mavx -I polybench/utilities polybench/utilities/polybench.c -I src/benchmark src/benchmark/bm.c
DUMP=-DPOLYBENCH_DUMP_ARRAYS
TIME=-DPOLYBENCH_TIME
EXTRA_FLAGS=-DPOLYBENCH_USE_C99_PROTO -D${SIZE}_DATASET ${TIME}

SRC_DIR=./src
BIN_DIR=./bin
TEST_DIR=./test
BINARY_DERICHE=${BIN_DIR}/deriche
BINARY_DERICHE_OMP=${BIN_DIR}/deriche_omp
BINARY_SEIDEL2D_OMP=${BIN_DIR}/seidel-2d_omp
BINARY_HEAT3D_OMP=${BIN_DIR}/heat-3d_omp
BINARY_HEAT3D_MPI=${BIN_DIR}/heat-3d_mpi
BINARY_HEAT3D_MPI_AVX2=${BIN_DIR}/heat-3d_mpi_avx2

OUTPUT_DERICHE=deriche_${SIZE}_out.txt

DERICHE_DIM=
ifdef DERICHE_DIM
	EXTRA_FLAGS := ${EXTRA_FLAGS} -DW=${DERICHE_DIM} -DH=${DERICHE_DIM}
endif

ifdef SEG_WIDTH
	EXTRA_FLAGS := ${EXTRA_FLAGS} -DSW=${SEG_WIDTH}
endif

ifdef DIM_W
ifdef DIM_H
		EXTRA_FLAGS := ${EXTRA_FLAGS} -DW=${DIM_W} -DH=${DIM_H}
endif
endif


$(shell mkdir -p bin)
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
CC=gcc
endif
ifeq ($(UNAME), Darwin)
CC=gcc-11
endif

MPI_CC=mpicc

all: deriche_omp deriche deriche_ref deriche_mpi_baseline deriche_mpi_rdma seidel-2d_omp heat-3d_omp heat-3d_mpi heat-3d_mpi_avx2

.PHONY: clean run

deriche_omp:
	${CC} ${CFLAGS} ${EXTRA_FLAGS} -fopenmp -o ${BINARY_DERICHE_OMP} ${SRC_DIR}/deriche/deriche_omp.c

deriche:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -Wno-unused-parameter -DNDEBUG -o ${BINARY_DERICHE} ${SRC_DIR}/deriche/deriche_mpi_r.c

deriche_ref:
	${CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/deriche_ref polybench/medley/deriche/deriche.c

deriche_mpi_baseline:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c

deriche_mpi_rdma:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c

deriche_mpi_segments:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c
	# for verification:
	# ${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c -DPOLYBENCH_DUMP_ARRAYS -DNO_BENCHMARK
	# mpiexec -np 2 ${BIN_DIR}/$@ > ${OUTPUT_DERICHE} 2>&1
	# diff ${OUTPUT_DERICHE} ${TEST_DIR}/deriche/deriche_${SIZE}_DATASET.out -s

deriche_mpi_avx:
	${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c
	# for verification:
	# ${MPI_CC} ${CFLAGS} ${EXTRA_FLAGS} -o ${BIN_DIR}/$@ ${SRC_DIR}/deriche/$@.c -DPOLYBENCH_DUMP_ARRAYS
	# mpiexec -np 2 ${BIN_DIR}/$@ > ${OUTPUT_DERICHE} 2>&1
	# diff ${OUTPUT_DERICHE} ${TEST_DIR}/deriche/deriche_${SIZE}_DATASET.out -s

run: deriche
	mpiexec -np 2 ./${BINARY_DERICHE} > ${OUTPUT_DERICHE} 2>&1
	diff ${OUTPUT_DERICHE} ${TEST_DIR}/deriche/deriche_${SIZE}_DATASET.out -s

seidel-2d_omp:
	${CC} ${CFLAGS}  -fopenmp -o ${BINARY_SEIDEL2D_OMP} ${SRC_DIR}/seidel-2d/seidel-2d_omp.c -I ${SRC_DIR}/seidel-2d ${EXTRA_FLAGS}

heat-3d_omp:
	${CC} ${CFLAGS}  -fopenmp -o ${BINARY_HEAT3D_OMP} ${SRC_DIR}/heat-3d/heat-3d_omp.c -I ${SRC_DIR}/heat-3d ${EXTRA_FLAGS}

heat-3d_mpi:
	${MPI_CC} ${CFLAGS}  ${SRC_DIR}/heat-3d/heat-3d_mpi.c  -o ${BINARY_HEAT3D_MPI} ${EXTRA_FLAGS}

heat-3d_mpi_avx2:
	${MPI_CC} ${CFLAGS}  ${SRC_DIR}/heat-3d/heat-3d_mpi_avx2.c  -o ${BINARY_HEAT3D_MPI_AVX2} ${EXTRA_FLAGS}

clean:
	@ rm -f ${BINARY_DERICHE} ${BINARY_DERICHE_OMP} ${BINARY_SEIDEL2D_OMP} ${BINARY_HEAT3D_OMP}
