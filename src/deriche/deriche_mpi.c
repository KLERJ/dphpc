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
#define DEBUG_PRINT_ARRAY(w, h, array) print_array(w, h, array);
#else
#define DEBUG_PRINT(format, args...)
#define DEBUG_PRINT_ARRAY(w, h, array)
#endif

#define ROOT_RANK 0

/* Include polybench common header. */
#include <polybench.h>

/* static void malloc_2d_array(int rows, int cols, DATA_TYPE **array) {} */

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

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int w, int h,
                        DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, imgOut[i][j]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}

static void deriche_horizontal(int w, int h, DATA_TYPE alpha,
                               DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
                               DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
                               DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
                               DATA_TYPE POLYBENCH_2D(y2, W, H, w, h)) {

  // Declare variables
  int i, j;
  DATA_TYPE xm1, ym1, ym2;
  DATA_TYPE xp1, xp2;
  /* DATA_TYPE tp1, tp2; */
  DATA_TYPE yp1, yp2;

  DATA_TYPE k;
  DATA_TYPE a1, a2, a3, a4;
  DATA_TYPE b1, b2, c1, c2;

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
      y1[i][j] = a1 * imgIn[i][j] + a2 * xm1 + b1 * ym1 + b2 * ym2;
      xm1 = imgIn[i][j];
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
      xp1 = imgIn[i][j];
      yp2 = yp1;
      yp1 = y2[i][j];
    }
  }

  // Intermediate image
  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      imgOut[i][j] = c1 * (y1[i][j] + y2[i][j]);
    }
  }

  return;
}
static void deriche_vertical(int w, int h, DATA_TYPE alpha,
                             DATA_TYPE POLYBENCH_2D(imgIn, W, H, w, h),
                             DATA_TYPE POLYBENCH_2D(imgOut, W, H, w, h),
                             DATA_TYPE POLYBENCH_2D(y1, W, H, w, h),
                             DATA_TYPE POLYBENCH_2D(y2, W, H, w, h)) {

  // Declare variables
  int i, j;
  DATA_TYPE tm1, ym1, ym2;
  /* DATA_TYPE xp1, xp2; */
  DATA_TYPE tp1, tp2;
  DATA_TYPE yp1, yp2;

  DATA_TYPE k;
  DATA_TYPE a5, a6, a7, a8;
  DATA_TYPE b1, b2, c2;

  /* #pragma scop */
  // Initialize variables
  k = (1.0 - exp(-alpha)) * (1.0 - exp(-alpha)) /
      (1.0 + 2.0 * alpha * exp(-alpha) - exp(2.0 * alpha));
  a5 = k;
  a6 = k * exp(-alpha) * (alpha - 1.0);
  a7 = k * exp(-alpha) * (alpha + 1.0);
  a8 = -k * exp(-2.0 * alpha);
  b1 = pow(2.0, -alpha);
  b2 = -exp(-2.0 * alpha);
  c2 = 1;

  // Vertical pass
  // NOTE: Fuse loops for better locality
  for (j = 0; j < h; j++) {
    tm1 = 0.0;
    ym1 = 0.0;
    ym2 = 0.0;
    for (i = 0; i < w; i++) {
      y1[i][j] = a5 * imgIn[i][j] + a6 * tm1 + b1 * ym1 + b2 * ym2;
      tm1 = imgIn[i][j];
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
      tp1 = imgIn[i][j];
      yp2 = yp1;
      yp1 = y2[i][j];
    }
  }

  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      imgOut[i][j] = c2 * (y1[i][j] + y2[i][j]);
    }
  }
}

int main(int argc, char **argv) {

  /* Initialize MPI and related variables */
  MPI_Init(&argc, &argv);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  DEBUG_PRINT("MPI initialized with size %d - rank %d\n", size, rank);

  /* Retrieve problem size. */
  int width = W;
  int height = H;

  // NOTE: We assume W >= H so we can reuse the y1 and y2 arrays
  assert(W >= H);

  DEBUG_PRINT("Started with w=%d h=%d\n", width, height);

  // Determine how many rows each process will receive
  // Each process receives up to n + 1 rows, where n = h // size
  // NOTE: Try to cache align rows
  int rows_per_processor = height / size;
  int mod_rows = height % size;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(imgIn, DATA_TYPE, W, H, width, height);
  POLYBENCH_2D_ARRAY_DECL(imgOut, DATA_TYPE, W, H, width, height);

  // local partition of arrays
  POLYBENCH_2D_ARRAY_DECL(imgInLocal, DATA_TYPE, W, H, width, height);
  POLYBENCH_2D_ARRAY_DECL(imgOutLocal, DATA_TYPE, W, H, width, height);
  POLYBENCH_2D_ARRAY_DECL(y1, DATA_TYPE, W, H, width, height);
  POLYBENCH_2D_ARRAY_DECL(y2, DATA_TYPE, W, H, width, height);

  /* Initialize input image */
  if (rank == ROOT_RANK) {
    DEBUG_PRINT("Rank %d Initializing input image\n", rank);
    init_array(width, height, &alpha, POLYBENCH_ARRAY(imgIn),
               POLYBENCH_ARRAY(imgOut));
    DEBUG_PRINT_ARRAY(width, height, POLYBENCH_ARRAY(imgIn));
    DEBUG_PRINT("Rank %d Done\n", rank);
  }

  /* TODO: Start timer. */
  /* polybench_start_instruments; */

  /* Run kernel. */

  /**** Partition matrix row-wise across ranks ****/

  // Number of elements to send to each processor
  int *sendcounts = malloc(sizeof(int) * size);
  // Displacement relative to sendbuf from which to send data to each processor
  int *displs = malloc(sizeof(int) * size);
  // Initialize displs and sendcounts
  int offset = 0;
  if (rank == ROOT_RANK) {
    DEBUG_PRINT("i  displs[i]   sendcounts[i]\n");
  }
  for (int i = 0; i < size; i++) {
    int n_rows = rows_per_processor;
    if (mod_rows > 0) {
      n_rows++;
      mod_rows--;
    }
    int elements = n_rows * width;
    sendcounts[i] = elements;
    displs[i] = offset;
    offset += elements;
    if (rank == ROOT_RANK) {
      DEBUG_PRINT("%i   %i    %i\n", i, displs[i], sendcounts[i]);
    }
  }

  // Scatter input image across ranks
  MPI_Scatterv(imgIn, sendcounts, displs, MPI_DOUBLE, imgInLocal,
               sendcounts[rank], MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);
  DEBUG_PRINT("Rank %d MPI_Scatterv\n", rank);

  // Determine how many rows we received
  // NOTE: Is there a point to determining this more accurately?
  int n_local_rows = sendcounts[rank] / width;
  DEBUG_PRINT("Rank %d got %d elements:\n", rank, sendcounts[rank]);
  DEBUG_PRINT_ARRAY(width, n_local_rows, POLYBENCH_ARRAY(imgInLocal));

  // Run horizontal pass
  deriche_horizontal(width, n_local_rows, alpha, POLYBENCH_ARRAY(imgInLocal),
                     POLYBENCH_ARRAY(imgOutLocal), POLYBENCH_ARRAY(y1),
                     POLYBENCH_ARRAY(y2));
  DEBUG_PRINT("Rank %d Horizontal pass\n", rank);

  // Gather outputs (imgOutLocal) from ranks to assemble the intermediate
  // output image (stored in imgIn)
  // recvcounts = sendcounts, displs stays the same, store in imgIn

  DEBUG_PRINT("Rank %d computed %d elements:\n", rank, sendcounts[rank]);
  DEBUG_PRINT_ARRAY(width, n_local_rows, POLYBENCH_ARRAY(imgOutLocal));

  MPI_Gatherv(imgOutLocal, sendcounts[rank], MPI_DOUBLE, imgIn, sendcounts,
              displs, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);
  DEBUG_PRINT("Rank %d MPI_Gatherv\n", rank);

  /* NOTE: We could transpose the matrix for better locality during the
   * pass computation. However, we will also need to transpose it again to
   * recover the correct result
   */

  /**** Partition matrix column-wise across ranks ****/

  // Determine how many columns each process will receive
  // Each process receives up to n + 1 columns, where n = w // size
  int cols_per_processor = width / size;
  int mod_cols = width % size;

  // Recompute send_counts and displs
  // NOTE: We send data at the granularity of *columns* now, so we need to
  // change the sendcounds accordingly
  offset = 0;
  for (int i = 0; i < size; i++) {
    int n_cols = cols_per_processor;
    if (mod_cols > 0) {
      n_cols++;
      mod_cols--;
    }
    int elements = n_cols;
    sendcounts[i] = elements;
    displs[i] = offset;
    offset += elements;
  }

  if (rank == ROOT_RANK) {
    DEBUG_PRINT_ARRAY(width, height, POLYBENCH_ARRAY(imgIn));
  }

  // This is a bit more involved than scattering row wise
  // New MPI Datatype to represent a column
  MPI_Datatype imgCol, imgColType, localCol, localColType;
  int n_rows = height;
  // The number of columns depends on the processor and whether we are sending
  // or receiving
  int n_local_cols = sendcounts[rank] / height;
  int n_cols = width;

  // Define the type for the *sending* processor
  if (rank == ROOT_RANK) {
    MPI_Type_vector(n_rows, 1, n_cols, MPI_DOUBLE, &imgCol);
    MPI_Type_commit(&imgCol);
    // Resize to make this work for columns (distance between elements)
    MPI_Type_create_resized(imgCol, 0, 1 * sizeof(DATA_TYPE), &imgColType);
    MPI_Type_commit(&imgColType);
  }

  // Define the type for the *receiving* processors
  MPI_Type_vector(n_rows, 1, n_local_cols, MPI_DOUBLE, &localCol);
  MPI_Type_commit(&localCol);
  // Resize to make this work for columns (distance between elements)
  MPI_Type_create_resized(localCol, 0, 1 * sizeof(DATA_TYPE), &localColType);
  MPI_Type_commit(&localColType);

  MPI_Scatterv(imgIn, sendcounts, displs, imgColType, imgInLocal,
               sendcounts[rank], localColType, ROOT_RANK, MPI_COMM_WORLD);

  DEBUG_PRINT("Rank %d MPI_Scatterv\n", rank);
  DEBUG_PRINT_ARRAY(n_cols, width, POLYBENCH_ARRAY(imgInLocal));
  // Run vertical pass
  deriche_vertical(n_local_cols, height, alpha, POLYBENCH_ARRAY(imgInLocal),
                   POLYBENCH_ARRAY(imgOutLocal), POLYBENCH_ARRAY(y1),
                   POLYBENCH_ARRAY(y2));
  DEBUG_PRINT("Rank %d Vertical pass\n", rank);
  // Gather outputs from ranks
  // recvcounts = sendcounts, displs stays the same, store in imgOut
  MPI_Gatherv(imgOutLocal, sendcounts[rank], localColType, imgOut, sendcounts,
              displs, imgColType, ROOT_RANK, MPI_COMM_WORLD);

  DEBUG_PRINT("Rank %d MPI_Gatherv\n", rank);

  // Free local arrays
  POLYBENCH_FREE_ARRAY(imgInLocal);
  POLYBENCH_FREE_ARRAY(imgOutLocal);
  POLYBENCH_FREE_ARRAY(y1);
  POLYBENCH_FREE_ARRAY(y2);
  free(sendcounts);
  free(displs);

  /* Stop and print timer. */
  /* polybench_stop_instruments; */
  /* polybench_print_instruments; */

  /* Print the output image to 1) check correctness and 2) prevent
   * dead-code elimination*/
  if (rank == ROOT_RANK) {
    DEBUG_PRINT("Rank %d Printing output\n", rank);
    print_array(width, height, POLYBENCH_ARRAY(imgOut));
  }

  /* Free input/output array */
  if (rank == ROOT_RANK) {
    POLYBENCH_FREE_ARRAY(imgIn);
    POLYBENCH_FREE_ARRAY(imgOut);
  }

  DEBUG_PRINT("MPI_Finalize\n");
  MPI_Finalize();

  return 0;
}
