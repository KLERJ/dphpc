/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* heat-3d.c: this file is part of PolyBench/C */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bm.h"
#include "heat-3d.h"
#include <immintrin.h>
#include <mpi.h>
#include <stdbool.h>

// #define DEBUG

#ifndef Z_DIM
#define Z_DIM 0
#endif
#ifndef Y_DIM
#define Y_DIM 0
#endif
#ifndef X_DIM
#define X_DIM 0
#endif

#define MPI_DATATYPE MPI_DOUBLE
#define IDX(ARRAY, X, Y, Z, NY, NZ)                                            \
  ((ARRAY) + ((X) * (NY) * (NZ)) + ((Y) * (NZ)) + (Z))
#define IDX2D(ARRAY, X, Y, NY) ((ARRAY) + ((X) * (NY)) + (Y))

#define UPDATE_STEP_GENERAL(SELF, LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK)       \
  (SCALAR_VAL(0.125) * (*(FRONT) + *(BACK) + *(TOP) + *(BOTTOM) + *(LEFT) +    \
                        *(RIGHT) + SCALAR_VAL(2.0) * *(SELF)))

#define UPDATE_STEP(A, I, J, K, NY, NZ)                                        \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ),                  \
      IDX(A, (I)-1, J, K, NY, NZ), IDX(A, I, (J) + 1, K, NY, NZ),              \
      IDX(A, I, (J)-1, K, NY, NZ), IDX(A, I, J, (K) + 1, NY, NZ),              \
      IDX(A, I, J, (K)-1, NY, NZ)))

#define UPDATE_STEP_LEFT_FACE(A, LEFT, I, J, K, NY, NZ)                        \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ),                  \
      IDX(A, (I)-1, J, K, NY, NZ), IDX(A, I, (J) + 1, K, NY, NZ), (LEFT),      \
      IDX(A, I, J, (K) + 1, NY, NZ), IDX(A, I, J, (K)-1, NY, NZ)))

#define UPDATE_STEP_RIGHT_FACE(A, RIGHT, I, J, K, NY, NZ)                      \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ),                  \
      IDX(A, (I)-1, J, K, NY, NZ), (RIGHT), IDX(A, I, (J)-1, K, NY, NZ),       \
      IDX(A, I, J, (K) + 1, NY, NZ), IDX(A, I, J, (K)-1, NY, NZ)))

#define UPDATE_STEP_TOP_FACE(A, TOP, I, J, K, NY, NZ)                          \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ),                  \
      IDX(A, (I)-1, J, K, NY, NZ), IDX(A, I, (J) + 1, K, NY, NZ),              \
      IDX(A, I, (J)-1, K, NY, NZ), IDX(A, I, J, (K) + 1, NY, NZ), (TOP)))

#define UPDATE_STEP_BOTTOM_FACE(A, BOTTOM, I, J, K, NY, NZ)                    \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ),                  \
      IDX(A, (I)-1, J, K, NY, NZ), IDX(A, I, (J) + 1, K, NY, NZ),              \
      IDX(A, I, (J)-1, K, NY, NZ), (BOTTOM), IDX(A, I, J, (K)-1, NY, NZ)))

#define UPDATE_STEP_FRONT_FACE(A, FRONT, I, J, K, NY, NZ)                      \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), IDX(A, (I) + 1, J, K, NY, NZ), (FRONT),         \
      IDX(A, I, (J) + 1, K, NY, NZ), IDX(A, I, (J)-1, K, NY, NZ),              \
      IDX(A, I, J, (K) + 1, NY, NZ), IDX(A, I, J, (K)-1, NY, NZ)))

#define UPDATE_STEP_BACK_FACE(A, BACK, I, J, K, NY, NZ)                        \
  (UPDATE_STEP_GENERAL(                                                        \
      IDX(A, I, J, K, NY, NZ), (BACK), IDX(A, (I)-1, J, K, NY, NZ),            \
      IDX(A, I, (J) + 1, K, NY, NZ), IDX(A, I, (J)-1, K, NY, NZ),              \
      IDX(A, I, J, (K) + 1, NY, NZ), IDX(A, I, J, (K)-1, NY, NZ)))

#define LEFT_NEIGHBOR 0
#define RIGHT_NEIGHBOR 1
#define TOP_NEIGHBOR 2
#define BOTTOM_NEIGHBOR 3
#define FRONT_NEIGHBOR 4
#define BACK_NEIGHBOR 5

/* Array initialization. */
static inline void init_array(int n, DATA_TYPE *A) {
  int i, j, k;
  // int c = 0;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        *IDX(A, i, j, k, n, n) = (DATA_TYPE)(i * i + j + (n - k)) * 10 / (n);
}

static inline void my_init_array(int nx, int ny, int nz, DATA_TYPE *A) {
  int i, j, k;
  // int c = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++) {
        *IDX(A, i, j, k, ny, nz) =
            (DATA_TYPE)(i * i + j + (nz - k)) * 10 / (nz);

        // *IDX(A, i, j, k, ny, nz) = (i + 1) * c++;
      }
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static inline void print_array(int n, DATA_TYPE *A)

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
        if ((i * n * n + j * n + k) % 20 == 0)
          fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, "%0.2lf ", *IDX(A, i, j, k, n, n));
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}

static inline void my_print_array(FILE *ptr, int nx, int ny, int nz,
                                  DATA_TYPE *A)

{
  int i, j, k;

  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++) {
        fprintf(ptr, "%6.2lf ", *IDX(A, i, j, k, ny, nz));
      }
      fprintf(ptr, "\n");
    }
    fprintf(ptr, "\n");
  }
}

static void my_local_init_array(int fullcube_sizes[3], int subcube_sizes[3],
                                int coords[3], DATA_TYPE *__restrict__ A) {
  // int nx = fullcube_sizes[0];
  // int ny = fullcube_sizes[1];
  int nz = fullcube_sizes[2];

  int sx = subcube_sizes[0];
  int sy = subcube_sizes[1];
  int sz = subcube_sizes[2];

  int px = coords[0];
  int py = coords[1];
  int pz = coords[2];

  int i, j, k, i_full, j_full, k_full;
  // int c = 0;
  for (i = 0; i < sx; i++)
    for (j = 0; j < sy; j++)
      for (k = 0; k < sz; k++) {
        i_full = i + px * sx;
        j_full = j + py * sy;
        k_full = k + pz * sz;

        *IDX(A, i, j, k, sy, sz) =
            (DATA_TYPE)(i_full * i_full + j_full + (nz - k_full)) * 10 / (nz);

        // *IDX(A, i, j, k, ny, nz) = (i + 1) * c++;
      }
}

#ifdef DEBUG

static void print_face(FILE *ptr, int ni, int nj, DATA_TYPE *A)

{
  int i, j;
  for (i = 0; i < ni; i++) {
    for (j = 0; j < nj; j++) {
      fprintf(ptr, "%6.2lf ", *IDX2D(A, i, j, nj));
    }
    fprintf(ptr, "\n");
  }
}

#endif

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static void compute_inner_step_kernel_heat_3d(int nx, int ny, int nz,
                                              DATA_TYPE *__restrict__ A,
                                              DATA_TYPE *__restrict__ B) {
  // Row major array
  // i loops over X
  // j loops over Y
  // k loops over Z (contiguous)

#if defined(OLD)

  // Constant arrays
  double a[4] = {0.125, 0.125, 0.125, 0.125};
  __m256d consts1 = _mm256_loadu_pd(a);

  double m1[4] = {0.125, 0.25, 0.125, 0};
  __m256d mask1 = _mm256_loadu_pd(m1);

  double m2[4] = {0, 0.125, 0.25, 0.125};
  __m256d mask2 = _mm256_loadu_pd(m2);

  int i, j, k;

  for (i = 2; i < nx - 2; i++) {
    for (j = 2; j < ny - 2; j++) {
      for (k = 2; k < nz - 6; k += 4) {
        // Note: there are comments below showing AVX register contents, from
        // highest bits to lowest.

        // Summing elements not from the z axis
        __m256d ip1 = _mm256_loadu_pd(IDX(A, i + 1, j, k, ny, nz));
        __m256d im1 = _mm256_loadu_pd(IDX(A, i - 1, j, k, ny, nz));
        __m256d jp1 = _mm256_loadu_pd(IDX(A, i, j + 1, k, ny, nz));
        __m256d jm1 = _mm256_loadu_pd(IDX(A, i, j - 1, k, ny, nz));

        __m256d kp1 =
            _mm256_loadu_pd(IDX(A, i, j, k + 1, ny, nz)); // k+4, k+3, k+2, k+1
        __m256d km1 =
            _mm256_loadu_pd(IDX(A, i, j, k - 1, ny, nz)); // k+2, k+1, k, k-1

        // The following can all be done in parallel
        __m256d r1 = _mm256_add_pd(ip1, im1);
        __m256d r2 = _mm256_add_pd(jp1, jm1);

        __m256d k1 = _mm256_mul_pd(km1, mask1);
        __m256d k2 = _mm256_mul_pd(km1, mask2);
        __m256d k3 = _mm256_mul_pd(kp1, mask1);
        __m256d k4 = _mm256_mul_pd(kp1, mask2);

        k1 = _mm256_hadd_pd(k1, k2);
        k3 = _mm256_hadd_pd(k3, k4);

        __m128d h1 = _mm256_extractf128_pd(k1, 1);
        __m128d h2 = _mm256_extractf128_pd(k3, 1);
        __m128d l1 = _mm256_extractf128_pd(k1, 0);
        __m128d l2 = _mm256_extractf128_pd(k3, 0);

        h1 = _mm_add_pd(h1, l1);
        h2 = _mm_add_pd(h2, l2);

        k1 = _mm256_set_m128d(h2, h1);

        r1 = _mm256_add_pd(r1, r2);
        r1 = _mm256_mul_pd(r1, consts1);

        r1 = _mm256_add_pd(r1, k1);

        _mm256_storeu_pd(IDX(B, i, j, k, ny, nz), r1);

        // Summing elements from the z axis
        /*double* q = (double*) &r1;
         *IDX(B, i, j, k, ny, nz) = q[0] + SCALAR_VAL(0.125)*(*IDX(A, i, j, k-1,
         *ny, nz) + *IDX(A, i, j, k+1, ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i,
         *j, k, ny, nz); IDX(B, i, j, k+1, ny, nz) = q[1] +
         *SCALAR_VAL(0.125)*(*IDX(A, i, j, k, ny, nz) + *IDX(A, i, j, k+2, ny,
         *nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+1, ny, nz); IDX(B, i, j,
         *k+2, ny, nz) = q[2] + SCALAR_VAL(0.125)*(*IDX(A, i, j, k+1, ny, nz) +
         **IDX(A, i, j, k+3, ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+2,
         *ny, nz);
         *IDX(B, i, j, k+3, ny, nz) = q[3] + SCALAR_VAL(0.125)*(*IDX(A, i, j,
         *k+2, ny, nz) + *IDX(A, i, j, k+4, ny, nz)) + SCALAR_VAL(0.25) *
         **IDX(A, i, j, k+3, ny, nz); */

        /**IDX(B, i, j, k, ny, nz) =  UPDATE_STEP(A, i, j, k, ny, nz);
         *IDX(B, i, j, k+1, ny, nz) =  UPDATE_STEP(A, i, j, k+1, ny, nz);
         *IDX(B, i, j, k+2, ny, nz) =  UPDATE_STEP(A, i, j, k+2, ny, nz);
         *IDX(B, i, j, k+3, ny, nz) =  UPDATE_STEP(A, i, j, k+3, ny, nz);*/
      }
      while (k < nz - 2) {
        *IDX(B, i, j, k, ny, nz) = UPDATE_STEP(A, i, j, k, ny, nz);
        k++;
      }
    }
  }

#elif defined(NEW)

  double a[4] = {0.125, 0.125, 0.125, 0.125};
  __m256d consts1 = _mm256_loadu_pd(a);

  double m1[4] = {0.125, 0.25, 0.125, 0};
  __m256d mask1 = _mm256_loadu_pd(m1);

  double m2[4] = {0, 0.125, 0.25, 0.125};
  __m256d mask2 = _mm256_loadu_pd(m2);

  int i, j, k;

  for (i = 1; i < nx - 1; i += 2) {
    for (j = 1; j < ny - 1; j += 2) {
      for (k = 1; k < nz - 5; k += 4) {
        // Note: there are comments below showing AVX register contents, from
        // highest bits to lowest.

        // Summing elements not from the z axis
        // __m256d ip1 =
        //     _mm256_loadu_pd(IDX(A, i + 1, j, k, ny, nz)); // Back column
        // __m256d im1 =
        //     _mm256_loadu_pd(IDX(A, i - 1, j, k, ny, nz)); // Front column
        // __m256d jp1 =
        //     _mm256_loadu_pd(IDX(A, i, j + 1, k, ny, nz)); // Right column
        // __m256d jm1 =
        //     _mm256_loadu_pd(IDX(A, i, j - 1, k, ny, nz)); // Left column

        /*
              b0 b1
           l1 c2 c3 r1
           l0 c0 c1 r0
              f0 f1


        */

        /********************
         * Compute C0
         * ******************/
        __m256d c0u =
            _mm256_loadu_pd(IDX(A, i, j, k - 1, ny, nz)); // Center upper
        __m256d c0l =
            _mm256_loadu_pd(IDX(A, i, j, k + 1, ny, nz)); // Center Lower

        __m256d c2u = _mm256_loadu_pd(IDX(A, i + 1, j, k - 1, ny, nz)); // Back
        __m256d c2l = _mm256_loadu_pd(IDX(A, i + 1, j, k + 1, ny, nz)); // Back
        __m256d c1u = _mm256_loadu_pd(IDX(A, i, j + 1, k - 1, ny, nz)); // Right
        __m256d c1l = _mm256_loadu_pd(IDX(A, i, j + 1, k + 1, ny, nz)); // Right

        // Convert centers to sides
        __m256d middle_upper = _mm256_permute4x64_pd(c1u, 0xf9); // Shift right
        __m256d middle_lower = _mm256_permute4x64_pd(c1l, 0x90); // Shift left
        __m256d c1middle =
            _mm256_blend_pd(middle_upper, middle_lower,
                            0x0c); // Blend the two registers together

        middle_upper = _mm256_permute4x64_pd(c2u, 0xf9); // Shift right
        middle_lower = _mm256_permute4x64_pd(c2l, 0x90); // Shift left
        __m256d c2middle =
            _mm256_blend_pd(middle_upper, middle_lower,
                            0x0c); // Blend the two registers together

        __m256d f0 = _mm256_loadu_pd(IDX(A, i - 1, j, k, ny, nz)); // Front
        __m256d l0 = _mm256_loadu_pd(IDX(A, i, j - 1, k, ny, nz)); // Left
        // Compute Sides
        __m256d total = _mm256_add_pd(f0, l0);
        total = _mm256_add_pd(total, c1middle);
        total = _mm256_add_pd(total, c2middle);
        total = _mm256_mul_pd(total, consts1);

        // Compute center
        __m256d k1 = _mm256_mul_pd(c0u, mask1);
        __m256d k2 = _mm256_mul_pd(c0u, mask2);

        __m256d sumU = _mm256_hadd_pd(k1, k2);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumU = _mm256_hadd_pd(sumU, sumU);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        __m256d k3 = _mm256_mul_pd(c0l, mask1);
        __m256d k4 = _mm256_mul_pd(c0l, mask2);

        __m256d sumL = _mm256_hadd_pd(k1, k2);
        sumL = _mm256_permute4x64_pd(
            sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumL = _mm256_hadd_pd(sumL, sumL);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        __m256d center = _mm256_blend_pd(sumU, sumL, 0x0c);
        total = _mm256_add_pd(total, center);

        // Todo Add center

        _mm256_storeu_pd(IDX(B, i, j, k, ny, nz), total);

        /********************
         * Compute C3
         * ******************/

        __m256d b1 = _mm256_loadu_pd(IDX(A, i + 2, j + 1, k, ny, nz)); // Back
        __m256d r1 = _mm256_loadu_pd(IDX(A, i + 1, j + 2, k, ny, nz)); // Right

        total = _mm256_add_pd(b1, r1);
        total = _mm256_add_pd(total, c1middle);
        total = _mm256_add_pd(total, c2middle);
        total = _mm256_mul_pd(total, consts1);

        __m256d c3u = _mm256_loadu_pd(
            IDX(A, i + 1, j + 1, k - 1, ny, nz)); // Center upper
        __m256d c3l = _mm256_loadu_pd(
            IDX(A, i + 1, j + 1, k + 1, ny, nz)); // Center Lower

        k1 = _mm256_mul_pd(c3u, mask1);
        k2 = _mm256_mul_pd(c3u, mask2);

        sumU = _mm256_hadd_pd(k1, k2);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumU = _mm256_hadd_pd(sumU, sumU);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        k3 = _mm256_mul_pd(c3l, mask1);
        k4 = _mm256_mul_pd(c3l, mask2);

        sumL = _mm256_hadd_pd(k1, k2);
        sumL = _mm256_permute4x64_pd(
            sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumL = _mm256_hadd_pd(sumL, sumL);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        center = _mm256_blend_pd(sumU, sumL, 0x0c);
        total = _mm256_add_pd(total, center);

        _mm256_storeu_pd(IDX(B, i + 1, j + 1, k, ny, nz), total);

        /********************
         * Compute C1
         * ******************/

        middle_upper = _mm256_permute4x64_pd(c0u, 0xf9); // Shift right
        middle_lower = _mm256_permute4x64_pd(c0l, 0x90); // Shift left
        __m256d c0middle =
            _mm256_blend_pd(middle_upper, middle_lower,
                            0x0c); // Blend the two registers together

        middle_upper = _mm256_permute4x64_pd(c3u, 0xf9); // Shift right
        middle_lower = _mm256_permute4x64_pd(c3l, 0x90); // Shift left
        __m256d c3middle =
            _mm256_blend_pd(middle_upper, middle_lower,
                            0x0c); // Blend the two registers together

        __m256d f1 = _mm256_loadu_pd(IDX(A, i - 1, j + 1, k, ny, nz)); // Front
        __m256d r0 = _mm256_loadu_pd(IDX(A, i, j + 2, k, ny, nz));     // Right

        total = _mm256_add_pd(f1, r0);
        total = _mm256_add_pd(total, c0middle);
        total = _mm256_add_pd(total, c3middle);
        total = _mm256_mul_pd(total, consts1);

        k1 = _mm256_mul_pd(c1u, mask1);
        k2 = _mm256_mul_pd(c1u, mask2);

        sumU = _mm256_hadd_pd(k1, k2);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumU = _mm256_hadd_pd(sumU, sumU);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        k3 = _mm256_mul_pd(c1l, mask1);
        k4 = _mm256_mul_pd(c1l, mask2);

        sumL = _mm256_hadd_pd(k1, k2);
        sumL = _mm256_permute4x64_pd(
            sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumL = _mm256_hadd_pd(sumL, sumL);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        center = _mm256_blend_pd(sumU, sumL, 0x0c);
        total = _mm256_add_pd(total, center);

        _mm256_storeu_pd(IDX(B, i, j + 1, k, ny, nz), total);

        /********************
         * Compute C2
         * ******************/

        __m256d b0 = _mm256_loadu_pd(IDX(A, i + 2, j, k, ny, nz));     // Back
        __m256d l1 = _mm256_loadu_pd(IDX(A, i + 1, j - 1, k, ny, nz)); // Left

        total = _mm256_add_pd(b0, l1);
        total = _mm256_add_pd(total, c0middle);
        total = _mm256_add_pd(total, c3middle);
        total = _mm256_mul_pd(total, consts1);

        k1 = _mm256_mul_pd(c2u, mask1);
        k2 = _mm256_mul_pd(c2u, mask2);

        sumU = _mm256_hadd_pd(k1, k2);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumU = _mm256_hadd_pd(sumU, sumU);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        k3 = _mm256_mul_pd(c2l, mask1);
        k4 = _mm256_mul_pd(c2l, mask2);

        sumL = _mm256_hadd_pd(k1, k2);
        sumL = _mm256_permute4x64_pd(
            sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
        sumL = _mm256_hadd_pd(sumL, sumL);
        sumU = _mm256_permute4x64_pd(
            sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

        center = _mm256_blend_pd(sumU, sumL, 0x0c);
        total = _mm256_add_pd(total, center);

        _mm256_storeu_pd(IDX(B, i + 1, j, k, ny, nz), total);
      }
      while (k < nz - 2) {
        *IDX(B, i, j, k, ny, nz) = UPDATE_STEP(A, i, j, k, ny, nz);
        k++;
      }
    }
  }
#else

#endif
}

static void compute_border_step_kernel_heat_3d(int nx, int ny, int nz,
                                               DATA_TYPE *__restrict__ A,
                                               DATA_TYPE *__restrict__ B,
                                               int neighbors[6],
                                               DATA_TYPE *neighborFaces[6]) {

  /*************************************************
   * Computing the border cases
   * **************************
   *
   * The funciton has the following if structure:
   *
   * front:
   *  top:
   *    FRONT-TOP-LEFT Corner
   *    FRONT-TOP Edge
   *    FRONT-TOP-RIGHT Corner
   *  FRONT Face
   *  bottom:
   *    FRONT-BOTTOM-LEFT Corner
   *    FRONT-BOTTOM Edge
   *    FRONT-BOTTOM-RIGHT Corner
   * back
   *  top:
   *    BACK-TOP-LEFT Corner
   *    BACK-TOP Edge
   *    BACK-TOP-RIGHT Corner
   *  BACK Face
   *
   *  bottom:
   *    BACK-BOTTOM-LEFT Corner
   *    BACK-BOTTOM Edge
   *    BACK-BOTTOM-RIGHT Corner
   *
   *
   * left
   *  FRONT-LEFT Edge
   *  LEFT Face
   *  BACK-LEFT Edge
   *
   *  LEFT-TOP Edge
   *  LEFT-BOTTOM Edge
   *
   * right
   *  FRONT-RIGHT Edge
   *  RIGHT Face
   *  BACK-RIGHT Edge
   *
   *  RIGHT-TOP Edge
   *  RIGHT-BOTTOM Edge
   *
   * top
   *  TOP FACE
   *
   * bottom
   *  BOTTOM FACE
   * ***********************************************/
  int x, y, z;

  // Work out which neighbors we have
  const bool left = neighbors[LEFT_NEIGHBOR] != MPI_PROC_NULL;
  const bool right = neighbors[RIGHT_NEIGHBOR] != MPI_PROC_NULL;
  const bool top = neighbors[TOP_NEIGHBOR] != MPI_PROC_NULL;
  const bool bottom = neighbors[BOTTOM_NEIGHBOR] != MPI_PROC_NULL;
  const bool front = neighbors[FRONT_NEIGHBOR] != MPI_PROC_NULL;
  const bool back = neighbors[BACK_NEIGHBOR] != MPI_PROC_NULL;

  const int lastx = nx - 1;
  const int lasty = ny - 1;
  const int lastz = nz - 1;

  /*************
   * FRONT
   * **********/
  if (front) {

    if (top) {

      // Compute FRONT-TOP-LEFT corner
      if (left) {
        *IDX(B, 0, 0, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, 0, 0, ny, nz), IDX(A, 1, 0, 0, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], 0, 0, nz),
            IDX(A, 0, 1, 0, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], 0, 0, nz),
            IDX(A, 0, 0, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], 0, 0, ny));
      }

      // Compute  FRONT-TOP Edge
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, 0, y, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, y, 0, ny, nz), IDX(A, 1, y, 0, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], y, 0, nz),
            IDX(A, 0, y + 1, 0, ny, nz), IDX(A, 0, y - 1, 0, ny, nz),
            IDX(A, 0, y, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], 0, y, ny));
      }

      // Compute FRONT-TOP-RIGHT corner
      if (right) {
        *IDX(B, 0, lasty, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, lasty, 0, ny, nz), IDX(A, 1, lasty, 0, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], lasty, 0, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], 0, 0, nz),
            IDX(A, 0, lasty - 1, 0, ny, nz), IDX(A, 0, lasty, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], 0, lasty, ny));
      }
    }

    // Compute FRONT Y-Z Face
    for (y = 1; y < ny - 1; y++) {
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, 0, y, z, ny, nz) = UPDATE_STEP_FRONT_FACE(
            A, IDX2D(neighborFaces[FRONT_NEIGHBOR], y, z, nz), 0, y, z, ny, nz);
      }
    }

    if (bottom) {

      // Compute FRONT-BOTTOM-LEFT corner
      if (left) {
        *IDX(B, 0, 0, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, 0, lastz, ny, nz), IDX(A, 1, 0, lastz, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], 0, lastz, nz),
            IDX(A, 0, 1, lastz, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], 0, lastz, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], 0, 0, ny),
            IDX(A, 0, 0, lastz - 1, ny, nz));
      }

      // Compute  FRONT-BOTTOM Edge
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, 0, y, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, y, lastz, ny, nz), IDX(A, 1, y, lastz, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], y, lastz, nz),
            IDX(A, 0, y + 1, lastz, ny, nz), IDX(A, 0, y - 1, lastz, ny, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], 0, y, ny),
            IDX(A, 0, y, lastz - 1, ny, nz));
      }

      // Compute FRONT-BOTTOM-RIGHT corner
      if (right) {
        *IDX(B, 0, lasty, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, lasty, lastz, ny, nz), IDX(A, 1, lasty, lastz, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], lasty, lastz, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], 0, lastz, nz),
            IDX(A, 0, lasty - 1, lastz, ny, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], 0, lasty, ny),
            IDX(A, 0, lasty, lastz - 1, ny, nz));
      }
    }
  }

  /*************
   * BACK
   * **********/
  if (back) {

    if (top) {
      // Compute BACK-TOP-LEFT corner
      if (left) {
        *IDX(B, lastx, 0, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, 0, 0, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], 0, 0, nz),
            IDX(A, lastx - 1, 0, 0, ny, nz), IDX(A, lastx, 1, 0, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], lastx, 0, nz),
            IDX(A, lastx, 0, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], lastx, 0, ny));
      }

      // Compute BACK-TOP Edges
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, lastx, y, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, y, 0, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], y, 0, nz),
            IDX(A, lastx - 1, y, 0, ny, nz), IDX(A, lastx, y + 1, 0, ny, nz),
            IDX(A, lastx, y - 1, 0, ny, nz), IDX(A, lastx, y, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], lastx, y, ny));
      }

      // Compute BACK-TOP-RIGHT corner
      if (right) {
        *IDX(B, lastx, lasty, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, lasty, 0, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], lasty, 0, nz),
            IDX(A, lastx - 1, lasty, 0, ny, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], lastx, 0, nz),
            IDX(A, lastx, lasty - 1, 0, ny, nz),
            IDX(A, lastx, lasty, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], lastx, lasty, ny));
      }
    }

    // Compute BACK Y-Z Face
    for (y = 1; y < ny - 1; y++) {
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, nx - 1, y, z, ny, nz) = UPDATE_STEP_BACK_FACE(
            A, IDX2D(neighborFaces[BACK_NEIGHBOR], y, z, nz), nx - 1, y, z, ny,
            nz);
      }
    }

    if (bottom) {
      // Compute BACK-BOTTOM-LEFT corner
      if (left) {
        *IDX(B, lastx, 0, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, 0, lastz, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], 0, lastz, nz),
            IDX(A, lastx - 1, 0, lastz, ny, nz),
            IDX(A, lastx, 1, lastz, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], lastx, lastz, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], lastx, 0, ny),
            IDX(A, lastx, 0, lastz - 1, ny, nz));
      }

      // Compute BACK-BOTTOM Edges
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, lastx, y, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, y, lastz, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], y, lastz, nz),
            IDX(A, lastx - 1, y, lastz, ny, nz),
            IDX(A, lastx, y + 1, lastz, ny, nz),
            IDX(A, lastx, y - 1, lastz, ny, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], lastx, y, ny),
            IDX(A, lastx, y, lastz - 1, ny, nz));
      }

      // Compute BACK-BOTTOM-RIGHT corner
      if (right) {
        *IDX(B, lastx, lasty, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, lasty, lastz, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], lasty, lastz, nz),
            IDX(A, lastx - 1, lasty, lastz, ny, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], lastx, lastz, nz),
            IDX(A, lastx, lasty - 1, lastz, ny, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], lastx, lasty, ny),
            IDX(A, lastx, lasty, lastz - 1, ny, nz));
      }
    }
  }

  /*************
   * LEFT
   * **********/
  if (left) {

    // Compute  FRONT-LEFT Edge
    if (front)
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, 0, 0, z, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, 0, z, ny, nz), IDX(A, 1, 0, z, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], 0, z, nz),
            IDX(A, 0, 1, z, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], 0, z, nz),
            IDX(A, 0, 0, z + 1, ny, nz), IDX(A, 0, 0, z - 1, ny, nz));
      }

    // Compute LEFT X-Z
    for (x = 1; x < nx - 1; x++) {
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, x, 0, z, ny, nz) = UPDATE_STEP_LEFT_FACE(
            A, IDX2D(neighborFaces[LEFT_NEIGHBOR], x, z, nz), x, 0, z, ny, nz);
      }
    }

    // Compute BACK-LEFT Edge
    if (back)
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, lastx, 0, z, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, 0, z, ny, nz), IDX(A, lastx - 1, 0, z, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], 0, z, nz),
            IDX(A, lastx, 1, z, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], lastx, z, nz),
            IDX(A, lastx, 0, z + 1, ny, nz), IDX(A, lastx, 0, z - 1, ny, nz));
      }

    // Compute LEFT-TOP Edge
    if (top)
      for (x = 1; x < nx - 1; x++) {
        *IDX(B, x, 0, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, x, 0, 0, ny, nz), IDX(A, x - 1, 0, 0, ny, nz),
            IDX(A, x + 1, 0, 0, ny, nz), IDX(A, x, 1, 0, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], x, 0, nz),
            IDX(A, x, 0, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], x, 0, ny));
      }

    // Compute LEFT-BOTTOM Edge
    if (bottom)
      for (x = 1; x < nx - 1; x++) {
        *IDX(B, x, 0, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, x, 0, lastz, ny, nz), IDX(A, x - 1, 0, lastz, ny, nz),
            IDX(A, x + 1, 0, lastz, ny, nz), IDX(A, x, 1, lastz, ny, nz),
            IDX2D(neighborFaces[LEFT_NEIGHBOR], x, lastz, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], x, 0, ny),
            IDX(A, x, 0, lastz - 1, ny, nz));
      }
  }

  /*************
   * RIGHT
   * **********/
  if (right) {

    // Compute  FRONT-RIGHT Edge
    if (front)
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, 0, lasty, z, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, 0, lasty, z, ny, nz), IDX(A, 1, lasty, z, ny, nz),
            IDX2D(neighborFaces[FRONT_NEIGHBOR], lasty, z, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], 0, z, nz),
            IDX(A, 0, lasty - 1, z, ny, nz), IDX(A, 0, lasty, z + 1, ny, nz),
            IDX(A, 0, lasty, z - 1, ny, nz));
      }

    // Compute RIGHT X-Z Face
    for (x = 1; x < nx - 1; x++) {
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, x, ny - 1, z, ny, nz) = UPDATE_STEP_RIGHT_FACE(
            A, IDX2D(neighborFaces[RIGHT_NEIGHBOR], x, z, nz), x, ny - 1, z, ny,
            nz);
      }
    }

    // Compute BACK-RIGHT Edge
    if (back)
      for (z = 1; z < nz - 1; z++) {
        *IDX(B, lastx, lasty, z, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, lastx, lasty, z, ny, nz),
            IDX2D(neighborFaces[BACK_NEIGHBOR], lasty, z, nz),
            IDX(A, lastx - 1, lasty, z, ny, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], lastx, z, nz),
            IDX(A, lastx, lasty - 1, z, ny, nz),
            IDX(A, lastx, lasty, z + 1, ny, nz),
            IDX(A, lastx, lasty, z - 1, ny, nz));
      }

    // Compute RIGHT-TOP Edge
    if (top)
      for (x = 1; x < nx - 1; x++) {
        *IDX(B, x, lasty, 0, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, x, lasty, 0, ny, nz), IDX(A, x - 1, lasty, 0, ny, nz),
            IDX(A, x + 1, lasty, 0, ny, nz), IDX(A, x, lasty - 1, 0, ny, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], x, 0, nz),
            IDX(A, x, lasty, 1, ny, nz),
            IDX2D(neighborFaces[TOP_NEIGHBOR], x, lasty, ny));
      }

    // Compute RIGHT-BOTTOM Edge
    if (bottom)
      for (x = 1; x < nx - 1; x++) {
        *IDX(B, x, lasty, lastz, ny, nz) = UPDATE_STEP_GENERAL(
            IDX(A, x, lasty, lastz, ny, nz),
            IDX(A, x - 1, lasty, lastz, ny, nz),
            IDX(A, x + 1, lasty, lastz, ny, nz),
            IDX(A, x, lasty - 1, lastz, ny, nz),
            IDX2D(neighborFaces[RIGHT_NEIGHBOR], x, lastz, nz),
            IDX2D(neighborFaces[BOTTOM_NEIGHBOR], x, lasty, ny),
            IDX(A, x, lasty, lastz - 1, ny, nz));
      }
  }

  /*************
   * TOP
   * **********/
  if (top) {
    // Compute TOP X-Y Face
    for (x = 1; x < nx - 1; x++) {
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, x, y, 0, ny, nz) = UPDATE_STEP_TOP_FACE(
            A, IDX2D(neighborFaces[TOP_NEIGHBOR], x, y, ny), x, y, 0, ny, nz);
      }
    }
  }

  /*************
   * BOTTOM
   * **********/
  if (bottom) {
    // Compute BOTTOM X-Y Face
    for (x = 1; x < nx - 1; x++) {
      for (y = 1; y < ny - 1; y++) {
        *IDX(B, x, y, nz - 1, ny, nz) = UPDATE_STEP_BOTTOM_FACE(
            A, IDX2D(neighborFaces[BOTTOM_NEIGHBOR], x, y, ny), x, y, nz - 1,
            ny, nz);
      }
    }
  }
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */

  int nx = N;
  int ny = N;
  int nz = N;
  int sx, sy, sz;
  int tsteps = TSTEPS;
  char *output_path;
  char *benchmark_path;

  bm_handle benchmark_exclusive_compute;
  bm_handle benchmark_overlap_compute;
  bm_handle benchmark_communication;
  bm_handle benchmark_iter;
  bm_handle benchmark_scatter_gather;
  bm_handle benchmark_total;

  if (argc == 3) {
    output_path = argv[1];
    benchmark_path = argv[2];
  } else if (argc == 5) {
    sscanf(argv[1], "%d", &nx);
    sscanf(argv[2], "%d", &tsteps);
    ny = nz = nx;
    output_path = argv[3];
    benchmark_path = argv[4];
  } else {
    fprintf(stderr,
            "Usage: %s <N_SIZE> <T_STEPS> OUTPUT_DUMP BENCHMARK_DUMP_PATH\n",
            argv[0]);
    exit(1);
  }

  bm_init(&benchmark_exclusive_compute, 2 * tsteps);
  bm_init(&benchmark_overlap_compute, 2 * tsteps);
  bm_init(&benchmark_communication, 2 * tsteps);
  bm_init(&benchmark_iter, 2 * tsteps);
  bm_init(&benchmark_scatter_gather, 2);
  bm_init(&benchmark_total, 1);

  DATA_TYPE *A, *B, *tmp;

  MPI_Init(&argc, &argv);
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // Number of processors along a dimension
  int pdims[3] = {0, 0, 0};
  MPI_Dims_create(p, 3, pdims);
  int px = pdims[0];
  int py = pdims[1];
  int pz = pdims[2];

#ifdef DEBUG

  char output_name[32];
  snprintf(output_name, 32, "output_%d.txt", rank);
  FILE *ptr = fopen(output_name, "w");

#endif

  if (rank == 0) {
    printf("Processors: %d N: %d, Steps:%d\n", p, nx, tsteps);
  }

  // Set up cubical topology topocomm
  int periods[3] = {X_DIM, Y_DIM, Z_DIM};
  MPI_Comm topocomm;
  MPI_Cart_create(comm, 3, pdims, periods, 0, &topocomm);

  // The current processor coordinates
  int coords[3];
  MPI_Cart_coords(topocomm, rank, 3, coords);

  printf("Rank: %d, x: %d, y: %d, z: %d\n", rank, coords[0], coords[1],
         coords[2]);

  // The sizes of the local volumes
  sx = nx / px;
  sy = ny / py;
  sz = nz / pz;

  int fullcube_sizes[3] = {nx, ny, nz};
  int subcube_sizes[3] = {sx, sy, sz};

  // TopBottom Z
  // LeftRight Y
  // FrontBack X
  int top, bottom, left, right, front, back;

  MPI_Cart_shift(topocomm, 0, 1, &front, &back);
  MPI_Cart_shift(topocomm, 1, 1, &left, &right);
  MPI_Cart_shift(topocomm, 2, 1, &top, &bottom);

  int neighbors[6] = {left, right, top, bottom, front, back};

  // Local cubes include face of neighbor
  A = (DATA_TYPE *)calloc(sx * sy * sz, sizeof(DATA_TYPE));
  B = (DATA_TYPE *)calloc(sx * sy * sz, sizeof(DATA_TYPE));

  if (A == NULL || B == NULL) {
    fprintf(stderr, "Error allocating A and B\n");
  }

  // Pointers to temporary neghbor faces
  DATA_TYPE *neighborFaces[6];

  for (int i = 0; i < 3; i++) {
    size_t size = 0;
    switch (i) {
    case 0:
      size = sx * sz; // Left  x-z
      break;          // Right x-z
    case 1:
      size = sx * sy; // Top x-y
      break;          // Bottom x-y
    case 2:
      size = sy * sz; // Front y-z
      break;          // Back y-z
    default:
      exit(EXIT_FAILURE);
    }

    neighborFaces[2 * i] = (DATA_TYPE *)calloc(size, sizeof(DATA_TYPE));
    neighborFaces[2 * i + 1] = (DATA_TYPE *)calloc(size, sizeof(DATA_TYPE));

    if (neighborFaces[2 * i] == NULL || neighborFaces[2 * i + 1] == NULL) {
      fprintf(stderr, "Error allocating neighbor arrays.\n");
      exit(EXIT_FAILURE);
    }
  }

#ifdef DEBUG

  fprintf(ptr, "rank: %d, px: %d, py: %d, pz: %d\n", rank, px, py, pz);
  fprintf(ptr, "rank: %d, rx: %d, ry: %d, rz: %d\n", rank, rx, ry, rz);
  fprintf(ptr, "rank: %d, sx: %d, sy: %d, sz: %d\n", rank, sx, sy, sz);

  fprintf(ptr, "    %4d%4d\n", top, back);
  fprintf(ptr, "%4d%4d%4d\n", left, rank, right);
  fprintf(ptr, "%4d%4d\n", front, bottom);
#endif

  MPI_Datatype subCubeSend, subCubeSend_tmp;
  MPI_Datatype subCubeRcv;

  MPI_Datatype xyFace;
  MPI_Datatype xzFace;
  MPI_Datatype yzFace;

  MPI_Datatype xyFace_neighbor;
  MPI_Datatype xzFace_neighbor;
  MPI_Datatype yzFace_neighbor;
  MPI_Datatype yRow;

  int starts[3] = {0, 0, 0};
  // Splitting Full Cube into subproblems
  MPI_Type_create_subarray(3, fullcube_sizes, subcube_sizes, starts,
                           MPI_ORDER_C, MPI_DATATYPE, &subCubeSend_tmp);
  MPI_Type_create_resized(subCubeSend_tmp, 0, sizeof(DATA_TYPE), &subCubeSend);

  // Subproblem receive cube
  MPI_Type_contiguous(sx * sy * sz, MPI_DATATYPE, &subCubeRcv);

  // XY-Face A[x][y][0] (double strided vectors)
  MPI_Type_vector(sy, 1, sz, MPI_DATATYPE, &yRow);
  MPI_Type_create_hvector(sx, 1, (sy * sz) * sizeof(MPI_DATATYPE), yRow,
                          &xyFace);

  // XZ-Face A[x][0][z] sx strided blocks of sz
  MPI_Type_vector(sx, sz, sz * sy, MPI_DATATYPE, &xzFace);

  // YZ-Face A[0][y][z] contiguous in sy*sz
  MPI_Type_contiguous(sy * sz, MPI_DATATYPE, &yzFace);

  // Contiguous neighbor Face buffers
  MPI_Type_contiguous(sx * sy, MPI_DATATYPE, &xyFace_neighbor);
  MPI_Type_contiguous(sx * sz, MPI_DATATYPE, &xzFace_neighbor);
  MPI_Type_contiguous(sy * sz, MPI_DATATYPE, &yzFace_neighbor);

  MPI_Type_commit(&subCubeSend);
  MPI_Type_commit(&subCubeRcv);

  MPI_Type_commit(&xyFace);
  MPI_Type_commit(&xzFace);
  MPI_Type_commit(&yzFace);

  MPI_Type_commit(&xyFace_neighbor);
  MPI_Type_commit(&xzFace_neighbor);
  MPI_Type_commit(&yzFace_neighbor);

  bm_start(&benchmark_scatter_gather);
#ifdef GLOBAL_INIT
  DATA_TYPE *full_cube = NULL;

  if (rank == 0) {

    /* Variable declaration/allocation. */
    full_cube = (DATA_TYPE *)calloc(nx * ny * nz, sizeof(DATA_TYPE));
    if (full_cube == NULL) {
      fprintf(stderr, "Falied allocating full cube memory.\n");
      exit(EXIT_FAILURE);
    }

    /* Initialize array(s). */
    // init_array(nx, full_cube);
    my_init_array(nx, ny, nz, full_cube);

    /* Start timer. */
    polybench_start_instruments;
  }

  /*************************************************************************
   *  Splittting problem into subcubes and scattering them among the procs
   *
   * ***********************************************************************/

  int *sendcounts = calloc(p, sizeof(int));
  int *displs = calloc(p, sizeof(int));
  if (sendcounts == NULL || displs == NULL) {
    fprintf(stderr, "Falied allocating sendcounts and displacement arrays\n");
    exit(EXIT_FAILURE);
  }
  int tmp_coords[3];
  for (int i = 0; i < p; i++) {
    sendcounts[i] = 1;

    MPI_Cart_coords(topocomm, i, 3, tmp_coords);
    displs[i] = ((tmp_coords[0] * sx) * (ny) * (nz)) +
                ((tmp_coords[1] * sy) * (nz)) + (tmp_coords[2] * sz);
  }

  MPI_Scatterv(full_cube, sendcounts, displs, subCubeSend, A, 1, subCubeRcv, 0,
               MPI_COMM_WORLD);
  memcpy(B, A, sizeof(DATA_TYPE) * sx * sy * sz);

#else

  my_local_init_array(fullcube_sizes, subcube_sizes, coords, A);

  memcpy(B, A, sizeof(DATA_TYPE) * sx * sy * sz);

  if (rank == 0) {
    polybench_start_instruments;
  }
#endif
  bm_stop(&benchmark_scatter_gather);

  /*************************************************************************
   *  Start timestep itteration
   *
   * ***********************************************************************/
  bm_start(&benchmark_total);
  for (int t = 1; t <= 2 * tsteps; t++) {

    bm_start(&benchmark_communication);
    bm_start(&benchmark_iter);

    // Send/Receive neighbor faces
    const int num_requsts = 2 * 6;
    MPI_Request requests[num_requsts];
    MPI_Isend(IDX(A, 0, 0, 0, sy, sz), 1, xzFace, left, 0, MPI_COMM_WORLD,
              requests + 0);
    MPI_Isend(IDX(A, 0, sy - 1, 0, sy, sz), 1, xzFace, right, 0, MPI_COMM_WORLD,
              requests + 1);

    MPI_Isend(IDX(A, 0, 0, 0, sy, sz), 1, xyFace, top, 0, MPI_COMM_WORLD,
              requests + 2);
    MPI_Isend(IDX(A, 0, 0, sz - 1, sy, sz), 1, xyFace, bottom, 0,
              MPI_COMM_WORLD, requests + 3);

    MPI_Isend(IDX(A, 0, 0, 0, sy, sz), 1, yzFace, front, 0, MPI_COMM_WORLD,
              requests + 4);
    MPI_Isend(IDX(A, sx - 1, 0, 0, sy, sz), 1, yzFace, back, 0, MPI_COMM_WORLD,
              requests + 5);

    MPI_Irecv(neighborFaces[0], 1, xzFace_neighbor, left, 0, MPI_COMM_WORLD,
              requests + 6);
    MPI_Irecv(neighborFaces[1], 1, xzFace_neighbor, right, 0, MPI_COMM_WORLD,
              requests + 7);

    MPI_Irecv(neighborFaces[2], 1, xyFace_neighbor, top, 0, MPI_COMM_WORLD,
              requests + 8);
    MPI_Irecv(neighborFaces[3], 1, xyFace_neighbor, bottom, 0, MPI_COMM_WORLD,
              requests + 9);

    MPI_Irecv(neighborFaces[4], 1, yzFace_neighbor, front, 0, MPI_COMM_WORLD,
              requests + 10);
    MPI_Irecv(neighborFaces[5], 1, yzFace_neighbor, back, 0, MPI_COMM_WORLD,
              requests + 11);

    bm_start(&benchmark_overlap_compute);
    // Start computing inner cube (which is independedt of the neighbors)
    compute_inner_step_kernel_heat_3d(sx, sy, sz, A, B);
    bm_stop(&benchmark_overlap_compute);

    // Wait for all the communication to complete
    MPI_Waitall(num_requsts, requests, MPI_STATUS_IGNORE);
    bm_stop(&benchmark_communication);

    bm_start(&benchmark_exclusive_compute);
    // Compute the border faces
    compute_border_step_kernel_heat_3d(sx, sy, sz, A, B, neighbors,
                                       neighborFaces);

    // Swap A and B
    tmp = A;
    A = B;
    B = tmp;

    bm_stop(&benchmark_exclusive_compute);
    bm_stop(&benchmark_iter);
    if (rank == 0 && t % 16 == 0) {
      printf("Iter %d\n", t);
    }
  }
  bm_stop(&benchmark_total);

  bm_start(&benchmark_scatter_gather);
#ifdef GLOBAL_INIT

  for (int i = 0; i < p; i++) {
    sendcounts[i] = 1;

    MPI_Cart_coords(topocomm, i, 3, tmp_coords);
    displs[i] = ((tmp_coords[0] * sx) * (ny) * (nz)) +
                ((tmp_coords[1] * sy) * (nz)) + (tmp_coords[2] * sz);
  }

  MPI_Gatherv(A, 1, subCubeRcv, full_cube, sendcounts, displs, subCubeSend, 0,
              MPI_COMM_WORLD);

  free(sendcounts);
  free(displs);
#else

#endif
  bm_stop(&benchmark_scatter_gather);

  /* Be clean. */
  if (rank == 0) {
    polybench_stop_instruments;
    polybench_print_instruments;

    /* Prevent dead-code elimination. All live-out data must be printed
       by the function call in argument. */
#ifdef GLOBAL_INIT

    FILE *out_ptr = fopen(output_path, "wb");

    fwrite(full_cube, sizeof(DATA_TYPE), nx * ny * nz, out_ptr);

    free(full_cube);
    fclose(out_ptr);
#else
    (void)output_path;
#endif
  }

#ifndef NO_BENCHMARKS // Only print benchmakrs, if we haven't disabled it
  // Dump benchmarks
  char bm_output_name[512];
  snprintf(bm_output_name, 512, "%s/benchmark_%d.csv", benchmark_path, rank);
  FILE *benchmark_dump = fopen(bm_output_name, "w");
  if (benchmark_dump == NULL) {
    fprintf(stderr, "Couldn't open file\n");
  }

  bm_print_events(&benchmark_total,
                  benchmark_dump); // total
  bm_print_events(&benchmark_scatter_gather,
                  benchmark_dump); // Scatter and gather benchmarks
  bm_print_events(&benchmark_communication,
                  benchmark_dump); // Communication time
  bm_print_events(&benchmark_overlap_compute,
                  benchmark_dump); // Overlap computation time
  bm_print_events(&benchmark_exclusive_compute,
                  benchmark_dump); // Exclusive compute (neighbors) time
  bm_print_events(&benchmark_iter, benchmark_dump); // Total time per iteration

  fclose(benchmark_dump);
#endif

  free(A);
  free(B);

#if DEBUG
  fclose(ptr);
#endif

  for (int i = 0; i < 6; i++) {
    free(neighborFaces[i]);
  }

  bm_destroy(&benchmark_communication);
  bm_destroy(&benchmark_overlap_compute);
  bm_destroy(&benchmark_exclusive_compute);
  bm_destroy(&benchmark_iter);

  MPI_Type_free(&subCubeSend);
  MPI_Type_free(&subCubeRcv);

  MPI_Type_free(&xyFace);
  MPI_Type_free(&xzFace);
  MPI_Type_free(&yzFace);

  MPI_Type_free(&xyFace_neighbor);
  MPI_Type_free(&xzFace_neighbor);
  MPI_Type_free(&yzFace_neighbor);

  MPI_Finalize();

  return 0;
}
