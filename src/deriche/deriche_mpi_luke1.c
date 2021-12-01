/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* deriche.c: this file is part of PolyBench/C */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "deriche.h"
#include <mpi.h>

/* Array initialization. */
static void init_array(int w, int h, DATA_TYPE *alpha,
                       DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
                       DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h)) {
  int i, j;

  *alpha = 0.25; // parameter of the filter

  // input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      imgIn[i][j] = (DATA_TYPE)((313 * i + 991 * j) % 65536) / 65535.0f;

}

static void my_init_array(int w, int h, DATA_TYPE *alpha, double* globalimgIn) {
  int i, j; 

  *alpha = 0.25;

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      globalimgIn[i*h + j] = (DATA_TYPE)((313 * i + 991 * j) % 65536) / 65535.0f;

}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int w, int h,
                        DATA_TYPE* imgOut)

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i*h + j]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Original code provided by Gael Deest */
static void kernel_deriche(int argc, char** argv, int w, int h, DATA_TYPE alpha) {
  DATA_TYPE k;
  DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
  DATA_TYPE b1, b2, c1, c2;

  alpha = 0.25; // TODO: Move the below to after this initialization in order to remove this
#pragma scop
  k = (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) *
      (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) /
      (SCALAR_VAL(1.0) + SCALAR_VAL(2.0) * alpha * EXP_FUN(-alpha) -
       EXP_FUN(SCALAR_VAL(2.0) * alpha));
  a1 = a5 = k;
  a2 = a6 = k * EXP_FUN(-alpha) * (alpha - SCALAR_VAL(1.0));
  a3 = a7 = k * EXP_FUN(-alpha) * (alpha + SCALAR_VAL(1.0));
  a4 = a8 = -k * EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  b1 = POW_FUN(SCALAR_VAL(2.0), -alpha);
  b2 = -EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  c1 = c2 = 1;

  MPI_Init(&argc, &argv);

  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Hi there. I am node %d of %d.\n", rank, nprocs);

  int i, j;
  DATA_TYPE xm1, tm1, ym1, ym2;
  DATA_TYPE xp1, xp2;
  DATA_TYPE tp1, tp2;
  DATA_TYPE yp1, yp2;

  // Arrays containing all input and output data.
  // Only for use by process 0.
  double* globalimgIn = NULL;
  double* globalimgOut = NULL;

  // Processor 0 initializes the global input image.
  if (rank == 0) {
    globalimgIn = (double*)malloc(w * h * sizeof(double));
    my_init_array(w, h, &alpha, globalimgIn);
    polybench_start_instruments;
  }

  double* localimgIn = (double*)malloc(h * w / nprocs * sizeof(double));
  double* y1 = (double*)malloc(h * w / nprocs * sizeof(double));
  double* y2 = (double*)malloc(h * w / nprocs * sizeof(double));
  double* localimgOuta = (double*)malloc(h * w / nprocs * sizeof(double));
  double* localimgOutb = (double*)malloc(h * w / nprocs * sizeof(double));

  // Processor 0 distributes rows of input image to workers.
  MPI_Scatter(globalimgIn, h * w / nprocs, MPI_DOUBLE, localimgIn, h * w / nprocs, MPI_DOUBLE,
          0, MPI_COMM_WORLD);

  for (i = 0; i < _PB_W/nprocs; i++) {
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    xm1 = SCALAR_VAL(0.0);
    for (j = 0; j < _PB_H; j++) {
      y1[i*_PB_H + j] = a1 *  localimgIn[i*_PB_H + j] + a2 * xm1 + b1 * ym1 + b2 * ym2;
      xm1 = localimgIn[i*_PB_H + j];
      ym2 = ym1;
      ym1 = y1[i*_PB_H + j];
    }
  }

  for (i = 0; i < _PB_W/nprocs; i++) {
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    xp1 = SCALAR_VAL(0.0);
    xp2 = SCALAR_VAL(0.0);
    for (j = _PB_H - 1; j >= 0; j--) {
      y2[i*h + j] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      xp2 = xp1;
      xp1 = localimgIn[i*h + j];
      yp2 = yp1;
      yp1 = y2[i*h + j];
    }
  }

  for (i = 0; i < _PB_W/nprocs; i++)
    for (j = 0; j < _PB_H; j++) {
      localimgOuta[i*h + j] = c1 * (y1[i*h + j] + y2[i*h + j]);
    }

  int idx = 0;
  for(int k = 0; k < nprocs; k++) {
      for (int i = 0; i < w/nprocs; i++) {
          for (int j = 0; j < h/nprocs; j++) {
              localimgOutb[idx] = localimgOuta[i*h + k*h/nprocs + j];
              idx++;
          }
      }
  }

  // Now need to exchange the localimgOut's between processors. 
  MPI_Alltoall(localimgOutb, h*w/nprocs/nprocs, MPI_DOUBLE, 
          localimgOuta, h*w/nprocs/nprocs, MPI_DOUBLE, MPI_COMM_WORLD);

  for (j = 0; j < _PB_H/nprocs; j++) {
    tm1 = SCALAR_VAL(0.0);
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    for (i = 0; i < _PB_W; i++) {
      y1[i*h/nprocs + j] = a5 * localimgOuta[i*h/nprocs + j] + a6 * tm1 + b1 * ym1 + b2 * ym2;
      tm1 = localimgOuta[i*h/nprocs + j];
      ym2 = ym1;
      ym1 = y1[i*h/nprocs + j];
    }
  }

  for (j = 0; j < _PB_H/nprocs; j++) {
    tp1 = SCALAR_VAL(0.0);
    tp2 = SCALAR_VAL(0.0);
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    for (i = _PB_W - 1; i >= 0; i--) {
      y2[i*h/nprocs + j] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;
      tp2 = tp1;
      tp1 = localimgOuta[i*h/nprocs + j];
      yp2 = yp1;
      yp1 = y2[i*h/nprocs + j];
    }
  }

  for (i = 0; i < _PB_W; i++) {
    for (j = 0; j < _PB_H/nprocs; j++) {
      localimgOuta[i*h/nprocs + j] = c2 * (y1[i*h/nprocs + j] + y2[i*h/nprocs + j]);
    }
  }

  if (rank == 0)
      globalimgOut = (double*)malloc(h * w * sizeof(double));

  MPI_Gather(localimgOuta, h*w/nprocs, MPI_DOUBLE, globalimgIn, h*w/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
      // Transposing
      idx = 0;
      for (int i = 0; i < w; i++) {
          for (int k = 0; k < nprocs; k++) {
              for (int j = 0; j < h/nprocs; j++) {
                  globalimgOut[idx] = globalimgIn[i*h/nprocs + k*w*h/nprocs + j];
                  idx++;
              }
          }
      }
      polybench_stop_instruments;
      polybench_print_instruments;
      polybench_prevent_dce(print_array(w, h, globalimgOut));
      free(globalimgIn);
      free(globalimgOut);
  }
  free(localimgIn);
  free(localimgOuta);
  free(localimgOutb);
  free(y1);
  free(y2);

  MPI_Finalize();

#pragma endscop
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  int w = W;
  int h = H;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;

  /* Run kernel. */
  kernel_deriche(argc, argv, w, h, alpha);


  return 0;
}
