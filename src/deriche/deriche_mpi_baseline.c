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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Include polybench common header. */
#include <mpi.h>
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bm.h"
#include "deriche.h"

#define ROOT_RANK 0
#define ind(i, j, size) ((i) * (size) + (j))

DATA_TYPE k;
DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
DATA_TYPE b1, b2, c1, c2;

// Benchmark handles
bm_handle benchmark_exclusive_compute;
bm_handle benchmark_overlap_compute;
bm_handle benchmark_communication;
bm_handle benchmark_iter;

/* Array initialization. */
static void init_array_private(long bw, long h, int rank, DATA_TYPE *imgIn) {
  long i, j;
  long istart = bw * rank;

  // input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < bw; i++)
    for (j = 0; j < h; j++)
      imgIn[ind(i, j, h)] =
          (DATA_TYPE)((313 * (istart + i) + 991 * j) % 65536) / 65535.0f;
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(long w, long h, DATA_TYPE *imgOut) {
  long i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER,
              imgOut[ind(i, j, h)]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}

static void deriche_horizontal(long bw, long h, DATA_TYPE *imgInPriv,
                               DATA_TYPE *y1) {
  DATA_TYPE xm1, ym1, ym2;
  DATA_TYPE xp1, xp2;
  DATA_TYPE yp1, yp2;
  for (long i = 0; i < bw; i++) {
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    xm1 = SCALAR_VAL(0.0);
    DATA_TYPE y1t;
    for (long j = 0; j < h; j++) {
      y1t = a1 * imgInPriv[ind(i, j, h)] + a2 * xm1 + b1 * ym1 + b2 * ym2;
      xm1 = imgInPriv[ind(i, j, h)];
      ym2 = ym1;
      ym1 = y1t;
      y1[ind(i, j, h)] = y1t;
    }
  }

  for (long i = 0; i < bw; i++) {
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    xp1 = SCALAR_VAL(0.0);
    xp2 = SCALAR_VAL(0.0);
    DATA_TYPE y2t;
    for (long j = h - 1; j >= 0; j--) {
      y2t = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      xp2 = xp1;
      xp1 = imgInPriv[ind(i, j, h)];
      yp2 = yp1;
      yp1 = y2t;
      y1[ind(i, j, h)] = c1 * (y1[ind(i, j, h)] + y2t);
    }
  }
}

static void deriche_vertical(long w, long bh, DATA_TYPE *imgOutPriv,
                             DATA_TYPE *y2) {

  DATA_TYPE tm1, ym1, ym2;
  DATA_TYPE tp1, tp2;
  DATA_TYPE yp1, yp2;

  for (long j = 0; j < bh; j++) {
    tm1 = SCALAR_VAL(0.0);
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    DATA_TYPE y1t;
    for (long i = 0; i < w; i++) {
      y1t = a5 * imgOutPriv[ind(i, j, bh)] + a6 * tm1 + b1 * ym1 + b2 * ym2;
      tm1 = imgOutPriv[ind(i, j, bh)];
      ym2 = ym1;
      ym1 = y1t;
      y2[ind(i, j, bh)] = y1t;
    }
  }

  for (long j = 0; j < bh; j++) {
    tp1 = SCALAR_VAL(0.0);
    tp2 = SCALAR_VAL(0.0);
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    DATA_TYPE y2t;
    for (long i = w - 1; i >= 0; i--) {
      y2t = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;
      tp2 = tp1;
      tp1 = imgOutPriv[ind(i, j, bh)];
      yp2 = yp1;
      yp1 = y2t;
      imgOutPriv[ind(i, j, bh)] = c2 * (y2[ind(i, j, bh)] + y2t);
    }
  }
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  long w = W;
  long h = H;

  // Initialize MPI environment.
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Let's assume size divides w and h for now.
  if (w % size != 0) {
    if (rank == ROOT_RANK)
      fprintf(stderr, "Image width %ld must be divisible by size %d.\n", w,
              size);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (h % size != 0) {
    if (rank == ROOT_RANK)
      fprintf(stderr, "Image height %ld must be divisible by size %d.\n", h,
              size);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  long bw = w / size;
  long bh = h / size;

  // Benchmark initialization
  char *benchmark_path = NULL;
  if (argc == 2) {
    benchmark_path = argv[1];
  }
  bm_init(&benchmark_exclusive_compute, 1);
  bm_init(&benchmark_overlap_compute, 1);
  bm_init(&benchmark_communication, 1);
  bm_init(&benchmark_iter, 1);

  MPI_Datatype _block_t;
  MPI_Type_vector(bw, bh, h, MPI_DOUBLE, &_block_t);
  MPI_Type_commit(&_block_t);
  MPI_Datatype block_t;
  MPI_Type_create_resized(_block_t, 0, bh * sizeof(double), &block_t);
  MPI_Type_commit(&block_t);

  MPI_Datatype _bh_cols_t;
  MPI_Type_vector(w, bh, h, MPI_DOUBLE, &_bh_cols_t);
  MPI_Type_commit(&_bh_cols_t);
  MPI_Datatype bh_cols_t;
  MPI_Type_create_resized(_bh_cols_t, 0, bh * sizeof(double), &bh_cols_t);
  MPI_Type_commit(&bh_cols_t);

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE *imgOut = NULL;
  DATA_TYPE *imgInPriv = malloc(bw * h * sizeof(DATA_TYPE));
  DATA_TYPE *imgOutPriv = malloc(w * bh * sizeof(DATA_TYPE));
  DATA_TYPE *y1 = malloc(bw * h * sizeof(DATA_TYPE));
  DATA_TYPE *y2 = malloc(w * bh * sizeof(DATA_TYPE));

  /* Initialize array(s). */
  if (rank == ROOT_RANK) {
    imgOut = malloc(w * h * sizeof(DATA_TYPE));
  }
  init_array_private(bw, h, rank, imgInPriv);

  alpha = 0.25; // parameter of the filter

  /* Start timer. */
  // polybench_start_instruments;
  double t_start = MPI_Wtime();
  bm_start(&benchmark_iter);
  bm_start(&benchmark_exclusive_compute);
  bm_start(&benchmark_overlap_compute);
  bm_stop(&benchmark_overlap_compute);

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

  /* Run kernel. */
  deriche_horizontal(bw, h, imgInPriv, y1);

  bm_pause(&benchmark_exclusive_compute);
  bm_start(&benchmark_communication);
  MPI_Alltoall(y1, 1, block_t, imgOutPriv, bw * bh, MPI_DOUBLE, MPI_COMM_WORLD);
  bm_stop(&benchmark_communication);
  bm_resume(&benchmark_exclusive_compute);

  deriche_vertical(w, bh, imgOutPriv, y2);

  /* Stop and print timer. */
  // polybench_stop_instruments;
  // polybench_print_instruments;
  double t_end = MPI_Wtime();
  bm_stop(&benchmark_exclusive_compute);
  bm_stop(&benchmark_iter);

  // MPI_Gather(imgOutPriv, w * bh, MPI_DOUBLE, imgOut, 1, bh_cols_t, ROOT_RANK,
  // MPI_COMM_WORLD);
  if (rank == ROOT_RANK) {
    printf("%0.6lf\n", t_end - t_start);
    /* Prevent dead-code elimination. All live-out data must be printed
       by the function call in argument. */
    /* polybench_prevent_dce(print_array(w, h, imgOut)); */
  }

  polybench_prevent_dce(print_array(w, bh, imgOutPriv));
  // Dump benchmarks
  if (benchmark_path != NULL) {
    char bm_output_name[512];
    snprintf(bm_output_name, 512, "%s/benchmark_%d.csv", benchmark_path, rank);
    FILE *benchmark_dump = fopen(bm_output_name, "w");
    if (benchmark_dump == NULL) {
      fprintf(stderr, "Couldn't open file\n");
    }

    bm_print_events(&benchmark_communication, benchmark_dump);
    bm_print_events(&benchmark_overlap_compute, benchmark_dump);
    bm_print_events(&benchmark_exclusive_compute, benchmark_dump);
    bm_print_events(&benchmark_iter, benchmark_dump);
    fclose(benchmark_dump);
  }

  /* Be clean. */
  bm_destroy(&benchmark_communication);
  bm_destroy(&benchmark_overlap_compute);
  bm_destroy(&benchmark_exclusive_compute);
  bm_destroy(&benchmark_iter);

  free(imgOut);
  free(imgInPriv);
  free(imgOutPriv);
  free(y1);
  free(y2);
  MPI_Type_free(&block_t);
  MPI_Type_free(&_block_t);
  MPI_Type_free(&bh_cols_t);
  MPI_Type_free(&_bh_cols_t);

  MPI_Finalize();
  return 0;
}
