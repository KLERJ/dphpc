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

#include <immintrin.h>

/* Include polybench common header. */
#include <mpi.h>
#include <polybench.h>

/* Include benchmark-specific header. */
#include "deriche.h"

#define ROOT_RANK 0
#define ind(i, j, size) ((i) * (size) + (j))

// How many rows does a segment contain?
#define SW 4

DATA_TYPE k;
__m256d k_;
DATA_TYPE a1, a2, a3, a4, a5, a6, a7, a8;
__m256d a1_, a2_, a3_, a4_, a5_, a6_, a7_, a8_;
DATA_TYPE b1, b2, c1, c2;
__m256d b1_, b2_, c1_, c2_;

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

static void deriche_horizontal(long bw, long h, long bh,
                               DATA_TYPE *restrict imgInPriv,
                               DATA_TYPE *restrict y1, int size,
                               MPI_Datatype segment_block_t) {

  DATA_TYPE ym1, ym2;
  for (long i = 0; i < bw; i+=1) {
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);

    DATA_TYPE y1t_0, y1t_1, y1t_2, y1t_3;

    // index 0 to 3
    long j = 0;
    long idx_0 = (i * h) + j;
    DATA_TYPE *firstHalf = (DATA_TYPE*) malloc(4*sizeof(DATA_TYPE));

    __m256d xm1_X1 = _mm256_load_pd(imgInPriv + idx_0);
    __m256d xm1_X2 = _mm256_set_pd(imgInPriv[idx_0 + 2], imgInPriv[idx_0 + 1], imgInPriv[idx_0], 0.0);
    __m256d a2xmX2 = _mm256_mul_pd(a2_, xm1_X2);
    __m256d firstHalf_ = _mm256_fmadd_pd(a1_, xm1_X1, a2xmX2);
    _mm256_store_pd(firstHalf, firstHalf_);

    y1t_0 = firstHalf[0] + b1 * ym1 + b2 * ym2;
    ym2 = ym1;
    ym1 = y1t_0;

    y1t_1 = firstHalf[1] + b1 * ym1 + b2 * ym2;
    ym2 = ym1;
    ym1 = y1t_1;

    y1t_2 = firstHalf[2] + b1 * ym1 + b2 * ym2;
    ym2 = ym1;
    ym1 = y1t_2;

    y1t_3 = firstHalf[3] + b1 * ym1 + b2 * ym2;
    ym2 = ym1;
    ym1 = y1t_3;

    y1[idx_0]     = y1t_0;
    y1[idx_0 + 1] = y1t_1;
    y1[idx_0 + 2] = y1t_2;
    y1[idx_0 + 3] = y1t_3;
    
    // index 4 to (h-1)
    for (long j = 4; j < h; j+=4) {
      long idx_0 = (i * h) + j;
    
      __m256d xm1_X1 = _mm256_load_pd(imgInPriv + idx_0);
      __m256d xm1_X2 = _mm256_load_pd(imgInPriv + idx_0 - 1);
      __m256d a2xmX2 = _mm256_mul_pd(a2_, xm1_X2);
      __m256d firstHalf_ = _mm256_fmadd_pd(a1_, xm1_X1, a2xmX2);
      _mm256_store_pd(firstHalf, firstHalf_);
     
      // xm1_1 = imgInPriv[idx_0];
      // y1t_0 = a1 * xm1_1 + a2 * xm1_0 + b1 * ym1_0 + b2 * ym2_0;
      // ym2_1 = ym1_0;
      // ym1_1 = y1t_0;

      // xm1_2 = imgInPriv[idx_1];
      // y1t_1 = a1 * xm1_2 + a2 * xm1_1 + b1 * ym1_1 + b2 * ym2_1;
      // ym2_2 = ym1_1;
      // ym1_2 = y1t_1;

      // xm1_3 = imgInPriv[idx_2];
      // y1t_2 = a1 * xm1_3 + a2 * xm1_2 + b1 * ym1_2 + b2 * ym2_2;
      // ym2_3 = ym1_2;
      // ym1_3 = y1t_2;

      // xm1_0 = imgInPriv[idx_3];
      // y1t_3 = a1 * xm1_0 + a2 * xm1_3 + b1 * ym1_3 + b2 * ym2_3;
      // ym2_0 = ym1_3;
      // ym1_0 = y1t_3;
      
      y1t_0 = firstHalf[0] + b1 * ym1 + b2 * ym2;
      ym2 = ym1;
      ym1 = y1t_0;

      y1t_1 = firstHalf[1] + b1 * ym1 + b2 * ym2;
      ym2 = ym1;
      ym1 = y1t_1;

      y1t_2 = firstHalf[2] + b1 * ym1 + b2 * ym2;
      ym2 = ym1;
      ym1 = y1t_2;

      y1t_3 = firstHalf[3] + b1 * ym1 + b2 * ym2;
      ym2 = ym1;
      ym1 = y1t_3;

      y1[idx_0]     = y1t_0;
      y1[idx_0 + 1] = y1t_1;
      y1[idx_0 + 2] = y1t_2;
      y1[idx_0 + 3] = y1t_3;
    }
  }

  DATA_TYPE yp1, yp2;
  for (long i = 0; i < bw; i+=1) {
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    DATA_TYPE y2t_0, y2t_1, y2t_2, y2t_3;

    // index h-1 to h-4
    long idx = (i * h) + (h - 1);
    DATA_TYPE *firstHalf = (DATA_TYPE*) malloc(4*sizeof(DATA_TYPE));

    __m256d xp1_X1 = _mm256_set_pd(imgInPriv[idx-2], imgInPriv[idx-1], imgInPriv[idx], 0.0);
    __m256d xp2_X2 = _mm256_set_pd(imgInPriv[idx-1], imgInPriv[idx], 0.0, 0.0);
    __m256d intrmd1 = _mm256_mul_pd(a4_, xp2_X2);
    __m256d firstHalf_ = _mm256_fmadd_pd(a3_, xp1_X1, intrmd1);
    _mm256_store_pd(firstHalf, firstHalf_);

    y2t_0 = firstHalf[0] + b1 * yp1 + b2 * yp2;
    yp2 = yp1;
    yp1 = y2t_0;

    y2t_1 = firstHalf[1] + b1 * yp1 + b2 * yp2;
    yp2 = yp1;
    yp1 = y2t_1;

    y2t_2 = firstHalf[2] + b1 * yp1 + b2 * yp2;
    yp2 = yp1;
    yp1 = y2t_2;

    y2t_3 = firstHalf[3] + b1 * yp1 + b2 * yp2;
    yp2 = yp1;
    yp1 = y2t_3;

    y1[idx]   = c1 * (y1[idx] + y2t_0);
    y1[idx-1] = c1 * (y1[idx-1] + y2t_1);
    y1[idx-2] = c1 * (y1[idx-2] + y2t_2);
    y1[idx-3] = c1 * (y1[idx-3] + y2t_3);

    // index h-5 to 0
    for (long j = h - 5; j >= 0; j-=4) {
      long idx = (i * h) + j;

      __m256d xp1_X1 = _mm256_set_pd(imgInPriv[idx-2], imgInPriv[idx-1], imgInPriv[idx], imgInPriv[idx+1]);
      __m256d xp2_X2 = _mm256_set_pd(imgInPriv[idx-1], imgInPriv[idx], imgInPriv[idx+1], imgInPriv[idx+2]);
      __m256d intrmd1 = _mm256_mul_pd(a4_, xp2_X2);
      __m256d firstHalf_ = _mm256_fmadd_pd(a3_, xp1_X1, intrmd1);
      _mm256_store_pd(firstHalf, firstHalf_);

      // y2t_0 = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      // xp2 = xp1;
      // xp1 = imgInPriv[idx];
      // yp2 = yp1;
      // yp1 = y2t_0;

      // y2t_1 = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      // xp2 = xp1;
      // xp1 = imgInPriv[idx-1];
      // yp2 = yp1;
      // yp1 = y2t_1;

      // y2t_2 = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      // xp2 = xp1;
      // xp1 = imgInPriv[idx-2];
      // yp2 = yp1;
      // yp1 = y2t_2;

      // y2t_3 = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      // xp2 = xp1;
      // xp1 = imgInPriv[idx-3];
      // yp2 = yp1;
      // yp1 = y2t_3;

      y2t_0 = firstHalf[0] + b1 * yp1 + b2 * yp2;
      yp2 = yp1;
      yp1 = y2t_0;

      y2t_1 = firstHalf[1] + b1 * yp1 + b2 * yp2;
      yp2 = yp1;
      yp1 = y2t_1;

      y2t_2 = firstHalf[2] + b1 * yp1 + b2 * yp2;
      yp2 = yp1;
      yp1 = y2t_2;

      y2t_3 = firstHalf[3] + b1 * yp1 + b2 * yp2;
      yp2 = yp1;
      yp1 = y2t_3;

      y1[idx]   = c1 * (y1[idx] + y2t_0);
      y1[idx-1] = c1 * (y1[idx-1] + y2t_1);
      y1[idx-2] = c1 * (y1[idx-2] + y2t_2);
      y1[idx-3] = c1 * (y1[idx-3] + y2t_3);
    }

    // Every SW rows (i.e. when one segment is complete) send it to processors
    if ((i != 0) && (i % SW == 0)) {
      for (int dst_rank = 0; dst_rank < size; dst_rank++) {
        MPI_Request req;
        size_t offset = (i - SW) * h + dst_rank * bh;
        MPI_Isend(y1 + offset, 1, segment_block_t, dst_rank, 0, MPI_COMM_WORLD,
                  &req);

        MPI_Request_free(&req);
        /* printf("Send to %d done offset: %zu\n", dst_rank, offset); */
      }
    }
  }

  // Send the last segment as well
  for (int dst_rank = 0; dst_rank < size; dst_rank++) {
    MPI_Request req;
    size_t offset = (bw - SW) * h + dst_rank * bh;
    MPI_Isend(y1 + offset, 1, segment_block_t, dst_rank, 0, MPI_COMM_WORLD,
              &req);

    MPI_Request_free(&req);
    /* printf("Send to %d done offset: %zu\n", dst_rank, offset); */
  }
}

static void deriche_vertical(long w, long bh, DATA_TYPE *restrict imgOutPriv,
                             DATA_TYPE *restrict y2,
                             MPI_Request *restrict requests) {

  DATA_TYPE tm1, ym1, ym2;
  DATA_TYPE tp1, tp2;
  DATA_TYPE yp1, yp2;

  for (long j = 0; j < bh; j++) {
    tm1 = SCALAR_VAL(0.0);
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    DATA_TYPE y1t;
    for (long i = 0; i < w; i++) {
      // Every SW rows, we need to make sure that the data we need is ready
      if (i % SW == 0) {
        MPI_Wait(requests + (i / SW), MPI_STATUS_IGNORE);
      }
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

  long bw = w / size; // == rows per processor == block width
  long bh = h / size; // == columns per processor == block height
  if (bw < SW) {
    fprintf(stderr, "Block width (%ld) less than segment width (%d)!\n", bw,
            SW);

    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  if (bw % SW != 0) {
    fprintf(stderr, "Block width (%ld) not divisible by segment width (%d)!\n",
            bw, SW);
    MPI_Abort(MPI_COMM_WORLD, 4);
  }

  // Segment block type : sw rows, bh columns each
  MPI_Datatype _segment_block_t, segment_block_t;
  MPI_Type_vector(SW, bh, h, MPI_DOUBLE, &_segment_block_t);
  MPI_Type_commit(&_segment_block_t);
  // Next segment block starts bh elements from the last one
  MPI_Type_create_resized(_segment_block_t, 0, bh * sizeof(double),
                          &segment_block_t);
  MPI_Type_commit(&segment_block_t);

  // Column segment type: w rows, bh columns each
  // Used only for gathering output
  MPI_Datatype _bh_cols_t, bh_cols_t;
  MPI_Type_vector(w, bh, h, MPI_DOUBLE, &_bh_cols_t);
  MPI_Type_commit(&_bh_cols_t);
  MPI_Type_create_resized(_bh_cols_t, 0, bh * sizeof(double), &bh_cols_t);
  MPI_Type_commit(&bh_cols_t);

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE *imgOut = NULL;
  DATA_TYPE *imgInPriv = malloc(bw * h * sizeof(DATA_TYPE));
  DATA_TYPE *y1 = malloc(bw * h * sizeof(DATA_TYPE));
  DATA_TYPE *y2 = malloc(w * bh * sizeof(DATA_TYPE));
  DATA_TYPE *imgOutPriv = malloc(w * bh * sizeof(DATA_TYPE));

  /* Initialize array(s). */
  if (rank == ROOT_RANK) {
    imgOut = malloc(w * h * sizeof(DATA_TYPE));
  }
  init_array_private(bw, h, rank, imgInPriv);

  alpha = 0.25; // parameter of the filter

  /* Start timer. */
  double t_start = MPI_Wtime();

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

  k_ = _mm256_set1_pd(k);
  a1_ = _mm256_set1_pd(a1);
  a2_ = _mm256_set1_pd(a2);
  a3_ = _mm256_set1_pd(a3);
  a4_ = _mm256_set1_pd(a4);
  a5_ = _mm256_set1_pd(a5);
  a6_ = _mm256_set1_pd(a6);
  a7_ = _mm256_set1_pd(a7);
  a8_ = _mm256_set1_pd(a8);
  b1_ = _mm256_set1_pd(b1);
  b2_ = _mm256_set1_pd(b2);
  c1_ = _mm256_set1_pd(c1);
  c2_ = _mm256_set1_pd(c2);

  /* Start receving segment blocks from other processors */
  const int segments_per_rank = bw / SW;
  const int n_requests = segments_per_rank * size; // == W / SW
  MPI_Request requests[n_requests];
  for (int src_rank = 0; src_rank < size; src_rank++) {
    // Receive segments from one rank
    for (int i = 0; i < segments_per_rank; i++) {
      size_t offset = src_rank * bw * bh + i * SW * bh;
      size_t req_offset = segments_per_rank * src_rank + i;
      size_t recv_count = SW * bh;
      /* printf("%d Irecv i: %d offset: %zu recvcount: %zu req_offset %zu\n",
       * rank, */
      /*        i, offset, recv_count, req_offset); */
      MPI_Irecv(imgOutPriv + offset, recv_count, MPI_DOUBLE, src_rank, 0,
                MPI_COMM_WORLD, requests + req_offset);
    }
  }

  /* Run kernel. */

  deriche_horizontal(bw, h, bh, imgInPriv, y1, size, segment_block_t);
  deriche_vertical(w, bh, imgOutPriv, y2, requests);

  /* Stop and print timer. */
  double t_end = MPI_Wtime();
#ifndef POLYBENCH_DUMP_ARRAYS
  printf("%d %0.6lf\n", rank, t_end - t_start);
#endif

#ifdef POLYBENCH_DUMP_ARRAYS
  MPI_Gather(imgOutPriv, w * bh, MPI_DOUBLE, imgOut, 1, bh_cols_t, ROOT_RANK,
             MPI_COMM_WORLD);
#endif
  if (rank == ROOT_RANK) {
    /* Prevent dead-code elimination. All live-out data must be printed
       by the function call in argument. */
    polybench_prevent_dce(print_array(w, h, imgOut));
  }

  /* Be clean. */
  free(imgOut);
  free(imgOutPriv);
  free(imgInPriv);
  free(y1);
  free(y2);
  MPI_Type_free(&segment_block_t);
  MPI_Type_free(&_segment_block_t);
  MPI_Type_free(&bh_cols_t);
  MPI_Type_free(&_bh_cols_t);

  MPI_Finalize();
  return 0;
}
