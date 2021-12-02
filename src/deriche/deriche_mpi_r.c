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
#define DEBUG_PRINT_ARRAY(w, h, array) debug_print_array(w, h, array);
#else
#define DEBUG_PRINT(format, args...)
#define DEBUG_PRINT_ARRAY(w, h, array)
#endif

#define IND(i, j, s) (i * s + j)

#define ROOT_RANK 0

int rank;
int size;

/* Include polybench common header. */
#include <polybench.h>

/* static void malloc_2d_array(int rows, int cols, DATA_TYPE **array) {} */

/* Array initialization. */
static void init_array(int w, int h, DATA_TYPE *imgIn) {
  int i, j;

  // input should be between 0 and 1 (grayscale image pixel)
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      imgIn[IND(i, j, w)] = (DATA_TYPE)((313 * i + 991 * j) % 65536) / 65535.0f;
  /* for debugging purposes: */
  /* imgIn[IND(i, j, h)] = (DATA_TYPE)(i * w) + j + 1; */
}

/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static void print_array(int w, int h, DATA_TYPE *imgOut)

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("imgOut");
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++) {
      if ((i * h + j) % 20 == 0)
        fprintf(POLYBENCH_DUMP_TARGET, "\n");
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER,
              imgOut[IND(i, j, h)]);
    }
  POLYBENCH_DUMP_END("imgOut");
  POLYBENCH_DUMP_FINISH;
}

// the same as print_array but prints row by row instead
static void debug_print_array(int w, int h, DATA_TYPE *imgOut)

{
  int i, j;

  fprintf(POLYBENCH_DUMP_TARGET, "debug_print ---");
  for (i = 0; i < h; i++) {
    fprintf(POLYBENCH_DUMP_TARGET, "\n %i:", i);
    for (j = 0; j < w; j++) {
      fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER,
              imgOut[IND(i, j, w)]);
    }
  }
  fprintf(POLYBENCH_DUMP_TARGET, "\n");
  fprintf(POLYBENCH_DUMP_TARGET, "debug_print ---\n");
}

static void deriche_horizontal(int w, int h, DATA_TYPE alpha, DATA_TYPE *imgIn,
                               DATA_TYPE *imgOut, DATA_TYPE *y1,
                               DATA_TYPE *y2) {

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
  k = (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) *
      (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) /
      (SCALAR_VAL(1.0) + SCALAR_VAL(2.0) * alpha * EXP_FUN(-alpha) -
       EXP_FUN(SCALAR_VAL(2.0) * alpha));
  a1 = k;
  a2 = k * EXP_FUN(-alpha) * (alpha - SCALAR_VAL(1.0));
  a3 = k * EXP_FUN(-alpha) * (alpha + SCALAR_VAL(1.0));
  a4 = -k * EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  b1 = POW_FUN(SCALAR_VAL(2.0), -alpha);
  b2 = -EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  c1 = c2 = 1;

  // Horizontal pass
  for (i = 0; i < w; i++) {
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    xm1 = SCALAR_VAL(0.0);
    for (j = 0; j < h; j++) {
      if (rank == 0) {
        DEBUG_PRINT("i=%i j=%i e=%i \n", i, j, IND(i, j, h));
      }
      y1[IND(i, j, h)] =
          a1 * imgIn[IND(i, j, h)] + a2 * xm1 + b1 * ym1 + b2 * ym2;
      xm1 = imgIn[IND(i, j, h)];
      ym2 = ym1;
      ym1 = y1[IND(i, j, h)];
    }
  }

  for (i = 0; i < w; i++) {
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    xp1 = SCALAR_VAL(0.0);
    xp2 = SCALAR_VAL(0.0);
    for (j = h - 1; j >= 0; j--) {
      y2[IND(i, j, h)] = a3 * xp1 + a4 * xp2 + b1 * yp1 + b2 * yp2;
      xp2 = xp1;
      xp1 = imgIn[IND(i, j, h)];
      yp2 = yp1;
      yp1 = y2[IND(i, j, h)];
    }
  }

  // Intermediate image
  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      imgOut[IND(i, j, h)] = c1 * (y1[IND(i, j, h)] + y2[IND(i, j, h)]);
    }
  }

  return;
}
static void deriche_vertical(int w, int h, DATA_TYPE alpha, DATA_TYPE *imgIn,
                             DATA_TYPE *imgOut, DATA_TYPE *y1, DATA_TYPE *y2) {

  // Declare variables
  int i, j;
  DATA_TYPE tm1, ym1, ym2;
  /* DATA_TYPE xp1, xp2; */
  DATA_TYPE tp1, tp2;
  DATA_TYPE yp1, yp2;

  DATA_TYPE k;
  DATA_TYPE a5, a6, a7, a8;
  DATA_TYPE b1, b2, c2;

  // Initialize variables
  k = (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) *
      (SCALAR_VAL(1.0) - EXP_FUN(-alpha)) /
      (SCALAR_VAL(1.0) + SCALAR_VAL(2.0) * alpha * EXP_FUN(-alpha) -
       EXP_FUN(SCALAR_VAL(2.0) * alpha));
  a5 = k;
  a6 = k * EXP_FUN(-alpha) * (alpha - SCALAR_VAL(1.0));
  a7 = k * EXP_FUN(-alpha) * (alpha + SCALAR_VAL(1.0));
  a8 = -k * EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  b1 = POW_FUN(SCALAR_VAL(2.0), -alpha);
  b2 = -EXP_FUN(SCALAR_VAL(-2.0) * alpha);
  c2 = 1;

  // Vertical pass
  for (j = 0; j < h; j++) {
    tm1 = SCALAR_VAL(0.0);
    ym1 = SCALAR_VAL(0.0);
    ym2 = SCALAR_VAL(0.0);
    for (i = 0; i < w; i++) {
      if (rank == 0) {
        DEBUG_PRINT("i=%i j=%i e=%i \n", i, j, IND(i, j, h));
      }
      y1[IND(i, j, h)] =
          a5 * imgIn[IND(i, j, h)] + a6 * tm1 + b1 * ym1 + b2 * ym2;
      tm1 = imgIn[IND(i, j, h)];
      ym2 = ym1;
      ym1 = y1[IND(i, j, h)];
    }
  }

  for (j = 0; j < h; j++) {
    tp1 = SCALAR_VAL(0.0);
    tp2 = SCALAR_VAL(0.0);
    yp1 = SCALAR_VAL(0.0);
    yp2 = SCALAR_VAL(0.0);
    for (i = w - 1; i >= 0; i--) {
      y2[IND(i, j, h)] = a7 * tp1 + a8 * tp2 + b1 * yp1 + b2 * yp2;
      tp2 = tp1;
      tp1 = imgIn[IND(i, j, h)];
      yp2 = yp1;
      yp1 = y2[IND(i, j, h)];
    }
  }

  for (i = 0; i < w; i++) {
    for (j = 0; j < h; j++) {
      imgOut[IND(i, j, h)] = c2 * (y1[IND(i, j, h)] + y2[IND(i, j, h)]);
    }
  }
}

int main(int argc, char **argv) {

  /* Initialize MPI and related variables */
  /* int size, rank; */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Retrieve problem size. */
  int width = W;
  int height = H;

  // NOTE: We assume W >= H so we can reuse the y1 and y2 arrays
  assert(W >= H);

  // Determine how many rows each process will receive

  // Each process receives up to n + 1 rows, where n = h // size
  // NOTE: Try to cache align rows
  int rows_per_processor = height / size;
  int mod_rows = height % size;

  // Number of elements to send to each processor
  int *sendcounts = malloc(sizeof(int) * size);
  // Displacement relative to sendbuf from which to send data to each processor
  int *displs = malloc(sizeof(int) * size);
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

  /* Variable declaration/allocation. */
  DATA_TYPE alpha = 0.25; // parameter of the filter
  DATA_TYPE *imgIn = NULL;
  DATA_TYPE *imgOut = NULL;
  int local_height = sendcounts[rank] / width;
  DATA_TYPE *imgInLocal = malloc(width * local_height * sizeof(DATA_TYPE));
  DATA_TYPE *imgOutLocal = malloc(width * local_height * sizeof(DATA_TYPE));
  DATA_TYPE *y1 = malloc(width * local_height * sizeof(DATA_TYPE));
  DATA_TYPE *y2 = malloc(width * local_height * sizeof(DATA_TYPE));

  // Allocate and initialize input/output on root
  if (rank == ROOT_RANK) {
    imgIn = malloc(width * height * sizeof(DATA_TYPE));
    imgOut = malloc(width * height * sizeof(DATA_TYPE));
    init_array(width, height, imgIn);
    DEBUG_PRINT_ARRAY(width, height, imgIn);
  }

  /* Start timer. */
  // TODO

  /* Run kernel. */

  /**** Partition matrix row-wise across ranks ****/

  // Scatter input image across ranks
  MPI_Scatterv(imgIn, sendcounts, displs, MPI_DOUBLE, imgInLocal,
               sendcounts[rank], MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

  // Run horizontal pass
  deriche_horizontal(width, local_height, alpha, imgInLocal, imgOutLocal, y1,
                     y2);

  // Gather outputs (imgOutLocal) from ranks to assemble the intermediate
  // output image (stored in imgOut)
  MPI_Gatherv(imgOutLocal, sendcounts[rank], MPI_DOUBLE, imgOut, sendcounts,
              displs, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);

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
    if (rank == ROOT_RANK) {
      DEBUG_PRINT("%i   %i    %i\n", i, displs[i], sendcounts[i]);
    }
  }

  /* if (rank == ROOT_RANK) { */
  /*   DEBUG_PRINT_ARRAY(width, height, imgIn); */
  /* } */
  if (rank == 0) {
    DEBUG_PRINT_ARRAY(width, local_height, imgInLocal);
  }

  // This is a bit more involved than scattering row wise
  // New MPI Datatype to represent a column
  MPI_Datatype imgCol, imgColType, localCol, localColType;
  // The number of columns depends on the processor and whether we are sending
  // or receiving
  int local_width = sendcounts[rank];

  // Define the type for the *sending* processor
  if (rank == ROOT_RANK) {
    MPI_Type_vector(height, 1, width, MPI_DOUBLE, &imgCol);
    MPI_Type_commit(&imgCol);
    // Resize to make this work for columns (distance between elements)
    MPI_Type_create_resized(imgCol, 0, 1 * sizeof(DATA_TYPE), &imgColType);
    MPI_Type_commit(&imgColType);
  }

  // Define the type for the *receiving* processors
  /* MPI_Type_vector(height, 1, n_local_cols, MPI_DOUBLE, &localCol); */
  MPI_Type_vector(height, 1, 1, MPI_DOUBLE, &localCol);
  MPI_Type_commit(&localCol);
  // Resize to make this work for columns (distance between elements)
  localColType = localCol;
  /* MPI_Type_create_resized(localCol, 0, height * sizeof(DATA_TYPE), */
  /* &localColType); */
  /* MPI_Type_commit(&localColType); */

  MPI_Scatterv(imgOut, sendcounts, displs, imgColType, imgInLocal,
               sendcounts[rank], localColType, ROOT_RANK, MPI_COMM_WORLD);

  if (rank == 0) {
    DEBUG_PRINT_ARRAY(height, local_width, imgInLocal);
  }
  // Run vertical pass
  deriche_horizontal(height, local_width, alpha, imgInLocal, imgOutLocal, y1,
                     y2);
  // Gather outputs from ranks
  // recvcounts = sendcounts, displs stays the same, store in imgOut

  MPI_Gatherv(imgOutLocal, sendcounts[rank], localColType, imgOut, sendcounts,
              displs, imgColType, ROOT_RANK, MPI_COMM_WORLD);

  /* Stop and print timer. */
  // TODO

  /* Print the output image to 1) check correctness and 2) prevent
   * dead-code elimination*/
  if (rank == ROOT_RANK) {
    DEBUG_PRINT_ARRAY(width, height, imgOut);
    print_array(width, height, imgOut);
  }

  /* Free arrays and types */
  free(imgInLocal);
  free(imgOutLocal);
  free(y1);
  free(y2);
  free(sendcounts);
  free(displs);
  MPI_Type_free(&localCol);
  if (rank == ROOT_RANK) {
    MPI_Type_free(&imgCol);
    MPI_Type_free(&imgColType);
    free(imgIn);
    free(imgOut);
  }

  MPI_Finalize();

  return 0;
}
