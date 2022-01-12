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
#include <omp.h>

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

/* Array initialization. */
// static void init_array(int n, DATA_TYPE POLYBENCH_3D(A, N, N, N, n, n, n),
//                        DATA_TYPE POLYBENCH_3D(B, N, N, N, n, n, n)) {
//   int i, j, k;

//   for (i = 0; i < n; i++)
//     for (j = 0; j < n; j++)
//       for (k = 0; k < n; k++)
//         A[i][j][k] = B[i][j][k] = (DATA_TYPE)(i + j + (n - k)) * 10 / (n);
// }

static inline void my_init_array(int nx, int ny, int nz, DATA_TYPE *A) {
  int i, j, k;
  int c = 0;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++) {
        // *IDX(A, i, j, k, ny, nz) =
        //(DATA_TYPE)(i * i + j + (nz - k)) * 10 / (nz);

        *IDX(A, i, j, k, ny, nz) = c++;
      }
}

// void print256_num(__m256d var) {
//   double val[4];
//   memcpy(val, &var, sizeof(val));
//   printf("%.2f %.2f %.2f %.2f\n", val[0], val[1], val[2], val[3]);
// }

// /* DCE code. Must scan the entire live-out data.
//    Can be used also to check the correctness of the output. */
// static void print_array(int n, DATA_TYPE POLYBENCH_3D(A, N, N, N, n, n, n))

// {
//   int i, j, k;

//   POLYBENCH_DUMP_START;
//   POLYBENCH_DUMP_BEGIN("A");
//   for (i = 0; i < n; i++)
//     for (j = 0; j < n; j++)
//       for (k = 0; k < n; k++) {
//         if ((i * n * n + j * n + k) % 20 == 0)
//           fprintf(POLYBENCH_DUMP_TARGET, "\n");
//         fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
//       }
//   POLYBENCH_DUMP_END("A");
//   POLYBENCH_DUMP_FINISH;
// }

bm_handle benchmark_exclusive_compute;
/* Main computational kernel. The whole function will be timed,
   including the call and return. */
void kernel_heat_3d(const int tsteps, const int n, DATA_TYPE *__restrict__ A,
                    DATA_TYPE *__restrict__ B) {

  const int nx = n, ny = n, nz = n;
  for (int t = 1; t <= 2 * tsteps; t++) {
    bm_start(&benchmark_exclusive_compute);
    //#pragma omp parallel for collapse(3)

#if defined(NEWNEW)

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

          __m256d kp1 = _mm256_loadu_pd(
              IDX(A, i, j, k + 1, ny, nz)); // k+4, k+3, k+2, k+1
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
           *IDX(B, i, j, k, ny, nz) = q[0] + SCALAR_VAL(0.125)*(*IDX(A, i, j,
           *k-1, ny, nz) + *IDX(A, i, j, k+1, ny, nz)) + SCALAR_VAL(0.25) *
           **IDX(A, i, j, k, ny, nz); IDX(B, i, j, k+1, ny, nz) = q[1] +
           *SCALAR_VAL(0.125)*(*IDX(A, i, j, k, ny, nz) + *IDX(A, i, j, k+2, ny,
           *nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+1, ny, nz); IDX(B, i, j,
           *k+2, ny, nz) = q[2] + SCALAR_VAL(0.125)*(*IDX(A, i, j, k+1, ny, nz)
           *+ *IDX(A, i, j, k+3, ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j,
           *k+2, ny, nz); IDX(B, i, j, k+3, ny, nz) = q[3] +
           *SCALAR_VAL(0.125)*(*IDX(A, i, j, k+2, ny, nz) + *IDX(A, i, j, k+4,
           *ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+3, ny, nz); */

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

#elif defined(OLD)

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

          __m256d kp1 = _mm256_loadu_pd(
              IDX(A, i, j, k + 1, ny, nz)); // k+4, k+3, k+2, k+1
          __m256d km1 =
              _mm256_loadu_pd(IDX(A, i, j, k - 1, ny, nz)); // k+2, k+1, k, k-1

          // The following can all be done in parallel
          __m256d r1 = _mm256_add_pd(ip1, im1);
          __m256d r2 = _mm256_add_pd(jp1, jm1);
          __m256d total = _mm256_add_pd(r1, r2);
          total = _mm256_mul_pd(total, consts1);

          __m256d k1 = _mm256_mul_pd(km1, mask1);
          __m256d k2 = _mm256_mul_pd(km1, mask2);
          __m256d k3 = _mm256_mul_pd(kp1, mask1);
          __m256d k4 = _mm256_mul_pd(kp1, mask2);

          __m256d sumU = _mm256_hadd_pd(k1, k2);

          sumU = _mm256_permute4x64_pd(
              sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumU = _mm256_hadd_pd(sumU, sumU);
          sumU = _mm256_permute4x64_pd(
              sumU, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

          __m256d sumL = _mm256_hadd_pd(k3, k4);

          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumL = _mm256_hadd_pd(sumL, sumL);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

          __m256d center = _mm256_blend_pd(sumU, sumL, 0x0c);

          total = _mm256_add_pd(total, center);

          _mm256_storeu_pd(IDX(B, i, j, k, ny, nz), total);

          // Summing elements from the z axis
          /*double* q = (double*) &r1;
           *IDX(B, i, j, k, ny, nz) = q[0] + SCALAR_VAL(0.125)*(*IDX(A, i, j,
           *k-1, ny, nz) + *IDX(A, i, j, k+1, ny, nz)) + SCALAR_VAL(0.25) *
           **IDX(A, i, j, k, ny, nz); IDX(B, i, j, k+1, ny, nz) = q[1] +
           *SCALAR_VAL(0.125)*(*IDX(A, i, j, k, ny, nz) + *IDX(A, i, j, k+2, ny,
           *nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+1, ny, nz); IDX(B, i, j,
           *k+2, ny, nz) = q[2] + SCALAR_VAL(0.125)*(*IDX(A, i, j, k+1, ny, nz)
           *+ *IDX(A, i, j, k+3, ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j,
           *k+2, ny, nz); IDX(B, i, j, k+3, ny, nz) = q[3] +
           *SCALAR_VAL(0.125)*(*IDX(A, i, j, k+2, ny, nz) + *IDX(A, i, j, k+4,
           *ny, nz)) + SCALAR_VAL(0.25) * *IDX(A, i, j, k+3, ny, nz); */

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

          __m256d c2u =
              _mm256_loadu_pd(IDX(A, i + 1, j, k - 1, ny, nz)); // Back
          __m256d c2l =
              _mm256_loadu_pd(IDX(A, i + 1, j, k + 1, ny, nz)); // Back
          __m256d c1u =
              _mm256_loadu_pd(IDX(A, i, j + 1, k - 1, ny, nz)); // Right
          __m256d c1l =
              _mm256_loadu_pd(IDX(A, i, j + 1, k + 1, ny, nz)); // Right

          // Convert centers to sides
          __m256d middle_upper =
              _mm256_permute4x64_pd(c1u, 0xf9); // Shift right
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

          __m256d sumL = _mm256_hadd_pd(k3, k4);

          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumL = _mm256_hadd_pd(sumL, sumL);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

          __m256d center = _mm256_blend_pd(sumU, sumL, 0x0c);

          total = _mm256_add_pd(total, center);

          _mm256_storeu_pd(IDX(B, i, j, k, ny, nz), total);

          /********************
           * Compute C3
           * ******************/

          __m256d b1 = _mm256_loadu_pd(IDX(A, i + 2, j + 1, k, ny, nz)); // Back
          __m256d r1 =
              _mm256_loadu_pd(IDX(A, i + 1, j + 2, k, ny, nz)); // Right

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

          sumL = _mm256_hadd_pd(k3, k4);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumL = _mm256_hadd_pd(sumL, sumL);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

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

          __m256d f1 =
              _mm256_loadu_pd(IDX(A, i - 1, j + 1, k, ny, nz));      // Front
          __m256d r0 = _mm256_loadu_pd(IDX(A, i, j + 2, k, ny, nz)); // Right

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

          sumL = _mm256_hadd_pd(k3, k4);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumL = _mm256_hadd_pd(sumL, sumL);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

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

          sumL = _mm256_hadd_pd(k3, k4);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2
          sumL = _mm256_hadd_pd(sumL, sumL);
          sumL = _mm256_permute4x64_pd(
              sumL, 0xd8); // Rearrange so that A1+B1, C1, B2, C2 + D2

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
#pragma omp parallel for collapse(3)
    for (int i = 1; i < nx - 1; i++) {
      for (int j = 1; j < ny - 1; j += 1) {
        for (int k = 1; k < nz - 1; k += 1) {
          *IDX(B, i, j, k, ny, nz) = UPDATE_STEP(A, i, j, k, ny, nz);
        }
      }
    }

#endif
    DATA_TYPE *tmp = B;
    B = A;
    A = tmp;
    bm_stop(&benchmark_exclusive_compute);
    if (t % 16 == 0) {
      printf("Iter %d\n", t);
    }
  }
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  int n;
  int tsteps;
  char *output_path;
  char *benchmark_path;

  if (argc == 3) {
    output_path = argv[1];
    benchmark_path = argv[2];
  } else if (argc == 5) {
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &tsteps);

    output_path = argv[3];
    (void)output_path;
    benchmark_path = argv[4];
  } else {
    fprintf(stderr,
            "Usage: %s <N_SIZE> <T_STEPS> OUTPUT_DUMP BENCHMARK_DUMP_PATH\n",
            argv[0]);
    exit(1);
  }
#pragma omp parallel
  {
    int i = omp_get_num_threads();
#pragma omp single
    { printf("OMP Num Threads: %d\n", i); }
  }

  /* Variable declaration/allocation. */
  DATA_TYPE *A = calloc(sizeof(DATA_TYPE), n * n * n);
  DATA_TYPE *B = calloc(sizeof(DATA_TYPE), n * n * n);

  /* Initialize array(s). */
  my_init_array(n, n, n, A);
  my_init_array(n, n, n, B);

  bm_init(&benchmark_exclusive_compute, 2 * tsteps);

#if defined(OLD)
  printf("Old AVX\n");
#elif defined(NEW)
  printf("New AVX\n");
#else
  printf("NO AVX\n");
#endif

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */

  kernel_heat_3d(tsteps, n, A, B);

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  char bm_output_name[1024];
  snprintf(bm_output_name, 1024, "%s/benchmark_0.csv", benchmark_path);
  FILE *benchmark_dump = fopen(bm_output_name, "w");
  if (benchmark_dump == NULL) {
    fprintf(stderr, "Couldn't open file\n");
  }

  bm_print_events(&benchmark_exclusive_compute,
                  benchmark_dump); // total
  bm_destroy(&benchmark_exclusive_compute);

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  // polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  free(A);
  free(B);

  return 0;
}
