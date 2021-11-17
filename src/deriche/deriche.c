#include "deriche.h"
#include <mpi.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef NDEBUG
#define DEBUG_PRINT(format, args...)                                           \
  printf("%s:%d %s() " format, __FILE__, __LINE__, __func__, ##args);
#elif
#define DEBUG_PRINT(format, args...)
#endif

#define ROOT_RANK 0

/* Include polybench common header. */
/* #include <polybench.h> */

/* static void malloc_2d_array(int rows, int cols, double **array) {} */

/* Array initialization. */
static void init_array(int w, int h, double *alpha, double input_image[][w]) {
  *alpha = 0.25; // parameter of the filter

  // input should be between 0 and 1 (grayscale image pixel)
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      input_image[i][j] = (double)((313 * i + 991 * j) % 65536) / 65535.0f;
    }
  }
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int width, int height, double **output_image) {
  int i, j;

  printf("output_image");
  for (i = 0; i < width; i++)
    for (j = 0; j < height; j++) {
      // Make it more readable by inserting newlines
      if ((i * height + j) % 20 == 0) {
        printf("\n");
      }
      printf("%lf ", output_image[i][j]);
    }
  printf("output image");
}

static void deriche_horizontal(int w, int h, double alpha,
                               double input_image[][w],
                               double output_image[][w], double y1[][w],
                               double y2[][w]) {

  // Declare variables
  int i, j;
  double xm1, tm1, ym1, ym2;
  double xp1, xp2;
  double tp1, tp2;
  double yp1, yp2;

  double k;
  double a1, a2, a3, a4;
  double b1, b2, c1, c2;

  // Initialize variables
  k = (1.0 - exp(-alpha)) * (1.0 - exp(-alpha)) /
      (1.0 + 2.0 * alpha * exp(-alpha) - exp(2.0 * alpha));
  a1 = k;
  a2 = k * exp(-alpha) * (alpha - 1.0);
  a3 = k * exp(-alpha) * (alpha + 1.0);
  a4 = -k * exp(-2.0 * alpha);
  b1 = pow(2.0, -alpha);
  b2 = -exp(-2.0 * alpha);
  c1 = c2 = 1;

  // Horizontal pass
  // NOTE: Fuse these loops for better locality
  for (i = 0; i < w; i++) {
    ym1 = 0.0;
    ym2 = 0.0;
    xm1 = 0.0;
    for (j = 0; j < h; j++) {
      y1[i][j] = a1 * input_image[i][j] + a2 * xm1 + b1 * ym1 + b2 * ym2;
      xm1 = input_image[i][j];
      ym2 = ym1;
      ym1 = y1[i][j];
    }
  }

  for (i = 0; i < w; i++) {
    yp1 = 0.0;
    yp2 = 0.0;
    xp1 = 0.0;
    xp2 = 0.0;
    for (j = h - 1; j >= 0; j--) {
      y2[i][j] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      xp2 = xp1;
      xp1 = input_image[i][j];
      yp2 = yp1;
      yp1 = y2[i][j];
    }
  }

  // Intermediate image
  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      output_image[i][j] = c1 * (y1[i][j] + y2[i][j]);
    }
  }

  return;
}

static void deriche_vertical(int w, int h, double alpha,
                             double input_image[][w], double output_image[][w],
                             double y1[][w], double y2[][w]) {

  // Declare variables
  int i, j;
  double xm1, tm1, ym1, ym2;
  double xp1, xp2;
  double tp1, tp2;
  double yp1, yp2;

  double k;
  double a5, a6, a7, a8;
  double b1, b2, c1, c2;

#pragma scop
  // Initialize variables
  k = (1.0 - exp(-alpha)) * (1.0 - exp(-alpha)) /
      (1.0 + 2.0 * alpha * exp(-alpha) - exp(2.0 * alpha));
  a5 = k;
  a6 = k * exp(-alpha) * (alpha - 1.0);
  a7 = k * exp(-alpha) * (alpha + 1.0);
  a8 = -k * exp(-2.0 * alpha);
  b1 = pow(2.0, -alpha);
  b2 = -exp(-2.0 * alpha);
  c1 = c2 = 1;

  // Vertical pass
  // NOTE: Fuse loops for better locality
  for (j = 0; j < h; j++) {
    tm1 = 0.0;
    ym1 = 0.0;
    ym2 = 0.0;
    for (i = 0; i < w; i++) {
      y1[i][j] = a5 * input_image[i][j] + a6 * tm1 + b1 * ym1 + b2 * ym2;
      tm1 = input_image[i][j];
      ym2 = ym1;
      ym1 = y1[i][j];
    }
  }

  for (j = 0; j < h; j++) {
    tp1 = 0.0;
    tp2 = 0.0;
    yp1 = 0.0;
    yp2 = 0.0;
    for (i = w - 1; i >= 0; i--) {
      y2[i][j] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;
      tp2 = tp1;
      tp1 = input_image[i][j];
      yp2 = yp1;
      yp1 = y2[i][j];
    }
  }

  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      output_image[i][j] = c2 * (y1[i][j] + y2[i][j]);
    }
  }
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  int width = W;
  int height = H;

  // NOTE: We assume W >= H so we can reuse the y1 and y2 arrays
  assert(W >= H);

  DEBUG_PRINT("Started with w=%d h=%d\n", width, height);

  /* Initialize MPI and related variables */
  int size;
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  DEBUG_PRINT("MPI initialized with size %d - rank %d\n", size, rank);

  // Determine how many rows each process will receive
  // Each process receives up to n + 1 rows, where n = h // size
  // NOTE: Try to cache align rows
  int rows_per_processor = height / size;
  int mod_rows = height % size;
  int elems_per_processor = (rows_per_processor + 1) * width;

  /* Variable declaration/allocation. */
  DEBUG_PRINT("Allocating arrays\n");
  double alpha;
  double(*input_image)[width] = malloc(sizeof(double[height][width]));
  double(*output_image)[width] = malloc(sizeof(double[height][width]));

  // local partition of arrays

  double(*input_image_local)[width] =
      malloc(sizeof(double[rows_per_processor + 1][width]));
  double(*output_image_local)[width] =
      malloc(sizeof(double[rows_per_processor + 1][width]));
  double(*y1)[width] = malloc(sizeof(double[rows_per_processor + 1][width]));
  double(*y2)[width] = malloc(sizeof(double[rows_per_processor + 1][width]));

  /* Initialize input image */
  DEBUG_PRINT("Initializing input image\n");
  init_array(width, height, &alpha, input_image);
  DEBUG_PRINT("Done\n");

  /* TODO: Start timer. */
  /* polybench_start_instruments; */

  /* Run kernel. */

  // Number of elements to send to each processor
  int sendcounts[size];
  // Displacement relative to sendbuf from which to send data to each processor
  int displs[size];
  // Initialize displs and sendcounts
  // TODO

  // Scatter input image across ranks
  DEBUG_PRINT("Rank %d MPI_Scatterv\n", rank);
  MPI_Scatterv(input_image, sendcounts, displs, MPI_DOUBLE, input_image_local,
               elems_per_processor, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

  // Determine how many rows we received
  // TODO
  int n_rows;

  // Run horizontal pass
  DEBUG_PRINT("Rank %d Horizontal pass\n", rank);
  deriche_horizontal(width, height, alpha, input_image_local,
                     output_image_local, y1, y2);

  // Gather outputs (output_image_local) from ranks to assemble the intermediate
  // output image (stored in input_image)
  DEBUG_PRINT("Rank %d MPI_Gatherv\n", rank);
  // recvcounts = sendcounts, displs stays the same, store in input_image
  MPI_Gatherv(output_image_local, elems_per_processor, MPI_DOUBLE, input_image,
              sendcounts, displs, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

  /* NOTE: We could transpose the matrix for better locality during the
   * pass computation. However, we will also need to transpose it again to
   * recover the correct result
   */

  // Determine how many columns each process will receive
  // Each process receives up to n + 1 columns, where n = w // size
  int cols_per_processor = width / size;
  int mod_cols = width % size;
  elems_per_processor = (cols_per_processor + 1) * width;

  // Recompute send_counts and displs

  // Scatter matrix column-wise across ranks
  DEBUG_PRINT("Rank %d MPI_Scatterv\n", rank);
  MPI_Scatterv(input_image, sendcounts, displs, MPI_DOUBLE, input_image_local,
               elems_per_processor, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

  // Run vertical pass
  DEBUG_PRINT("Rank %d Vertical pass\n", rank);
  deriche_vertical(width, height, alpha, input_image_local, output_image_local,
                   y1, y2);

  // Gather outputs from ranks
  DEBUG_PRINT("Rank %d MPI_Gatherv\n", rank);
  // recvcounts = sendcounts, displs stays the same, store in output_image
  MPI_Gatherv(output_image_local, elems_per_processor, MPI_DOUBLE, output_image,
              sendcounts, displs, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

  // Free local arrays
  free(input_image_local);
  free(output_image_local);
  free(input_image);
  free(y1);
  free(y2);
  /* Stop and print timer. */
  /* polybench_stop_instruments; */
  /* polybench_print_instruments; */

  /* Print the output image to 1) check correctness and 2) prevent
   * dead-code elimination*/

  /* Free output array */
  free(output_image);

  DEBUG_PRINT("MPI_Finalize\n");
  MPI_Finalize();

  return 0;
}
