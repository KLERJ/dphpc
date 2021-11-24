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

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "heat-3d.h"
#include <mpi.h>


#define MPI_DATATYPE MPI_DOUBLE
#define IDX(ARRAY, X, Y ,Z, NY, NZ) ((ARRAY) + ((X) * (NY) * (NZ)) + ((Y) * (NZ)) + (Z))



/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE *A)
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        *IDX(A, i, j, k, n, n) = (DATA_TYPE) (i + j*j + (n-k))* 10 / (n);
}

static
void my_init_array (int nx, int ny, int nz,
		 DATA_TYPE *A)
{
  int i, j, k;
  int c = 0;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++){
        // c = rand()%100;
        *IDX(A, i, j, k, ny, nz) =  (DATA_TYPE) (i + j*j + (nz-k))* 10 / (nz);
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE *A)

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
         if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(POLYBENCH_DUMP_TARGET, "%0.2lf ", *IDX(A, i, j, k, n , n));
      }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}

static
void my_print_array(int nx, int ny, int nz,
		 DATA_TYPE *A)

{
  int i, j, k;

  // POLYBENCH_DUMP_START;
  // POLYBENCH_DUMP_BEGIN("A");
  
  
  for (i = 0; i < nx; i++) {
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++) {
        //  if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
         fprintf(stderr, DATA_PRINTF_MODIFIER, *IDX(A, i, j, k, ny, nz));
      }
      fprintf(stderr, "\n");
    }
      fprintf(stderr, "\n\n");
  }
  // POLYBENCH_DUMP_END("A");
  // POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void compute_step_kernel_heat_3d(int nx,
          int ny,
          int nz,
		      DATA_TYPE *A,
		      DATA_TYPE *B)
{
  // Row major array
  // i loops over X
  // j loops over Y 
  // k loops over Z (contiguous)

  for (int i = 1; i < nx-1; i++) {
      for (int j = 1; j < ny-1; j++) {
          for (int k = 1; k < nz-1; k++) {    
              *IDX(B, i, j, k, ny, nz) =  SCALAR_VAL(0.125) * (*IDX(A, i+1, j, k, ny, nz) - SCALAR_VAL(2.0) * *IDX(A, i, j, k, ny, nz) + *IDX(A, i-1, j, k, ny, nz))
                                        + SCALAR_VAL(0.125) * (*IDX(A, i, j+1, k, ny, nz) - SCALAR_VAL(2.0) * *IDX(A, i, j, k, ny, nz) + *IDX(A, i, j-1, k, ny, nz))
                                        + SCALAR_VAL(0.125) * (*IDX(A, i, j, k+1, ny, nz) - SCALAR_VAL(2.0) * *IDX(A, i, j, k, ny, nz) + *IDX(A, i, j, k-1, ny, nz))
                                        + *IDX(A, i, j, k, ny, nz);
          }
      }
  }
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int nx = N;
  int ny = N;
  int nz = N;
  int tsteps = TSTEPS;

  (void)&print_array;
  (void)&init_array;

  double *full_cube = NULL, *A, *B;
 
  MPI_Init(&argc, &argv); 
  int rank, p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &p);

  // Number of processors along a dimension
  int pdims[3]={0,0,0};
  MPI_Dims_create(p, 3, pdims);
  int px = pdims[0];
  int py = pdims[1];
  int pz = pdims[2];

  if(rank == 0){
    fprintf(stdout, "Num: %d Procs: %d, Steps:%d\n", nx, p, tsteps);
  }

  printf("rank: %d, px: %d, py: %d, pz: %d\n", rank, px, py, pz);

  int periods[3] = {0, 0, 0};
  MPI_Comm topocomm;
  MPI_Cart_create(comm, 3, pdims, periods, 0, &topocomm);

  int coords[3];
  MPI_Cart_coords(topocomm, rank, 3, coords);

  // The current processor coordinates
  int rx = coords[0];
  int ry = coords[1];
  int rz = coords[2];

  printf("rank: %d, rx: %d, ry: %d, rz: %d\n", rank, rx, ry, rz);

  // The sizes of the local volumes
  int sizes[3];
  sizes[0] = nx / px;
  sizes[1] = ny / py;
  sizes[2] = nz / pz;

  
  int sx = sizes[0];
  int sy = sizes[1];
  int sz = sizes[2];
  printf("rank: %d, rxsx: %d, rysy: %d, rzsz: %d\n", rank, rx*sx, ry*sy, rz*sz);

  // Sizes of the local volumes including neighbor boundary
  int ssx = sx + 2;
  int ssy = sy + 2;
  int ssz = sz + 2;
  

  printf("rank: %d, sx: %d, sy: %d, sz: %d\n", rank, sx, sy, sz);
  
  // TopBottom Z
  // LeftRight Y
  // FrontBack X
  int top, bottom, left, right, front, back;

  MPI_Cart_shift(topocomm, 0, 1, &front, &back);
  MPI_Cart_shift(topocomm, 1, 1, &left, &right);
  MPI_Cart_shift(topocomm, 2, 1, &top, &bottom);

  
  // Local cubes include face of neighbor
  A = (double*)calloc(ssx* ssy * ssz, sizeof(double));
  B = (double*)calloc(ssx* ssy * ssz, sizeof(double));

  int orig_sizes[3] = {N, N, N};
  MPI_Datatype subCubeSend, subCubeSend_tmp;
  MPI_Datatype subCubeRcv;

  MPI_Datatype xyFace;
  MPI_Datatype xzFace;
  MPI_Datatype yzFace;
  MPI_Datatype yRow;


  int starts[3] = {0,0,0};
  MPI_Type_create_subarray(3, orig_sizes, sizes, starts, MPI_ORDER_C, MPI_DATATYPE, &subCubeSend_tmp);
  MPI_Type_create_resized(subCubeSend_tmp, 0,  sizeof(double),  &subCubeSend);

  MPI_Datatype y_s_slice;
  MPI_Type_vector(sy, sz, ssz, MPI_DATATYPE, &y_s_slice);
  MPI_Type_create_hvector(sx, 1, (ssy * ssz) * sizeof(MPI_DATATYPE), y_s_slice, &subCubeRcv);

  // A[x][y][0] strided with sx*sy elements, each of length 1. Offset sz  
  MPI_Type_vector(sy, 1, ssz, MPI_DATATYPE, &yRow);
  MPI_Type_create_hvector(sx, 1, (ssy * ssz) * sizeof(MPI_DATATYPE), yRow, &xyFace);

  // A[x][0][z] strided with ssx*ssz elements in blocks of sz, strided by ssz*ssy
  MPI_Type_vector(sx, sz, ssz*ssy, MPI_DATATYPE, &xzFace);
  
  // A[0][y][z] sy* sz elements in blocks of sz strided by ssz
  MPI_Type_vector(sy, sz, ssz,  MPI_DATATYPE, &yzFace);
  

  
  MPI_Type_commit(&subCubeSend);
  MPI_Type_commit(&subCubeRcv);

  MPI_Type_commit(&xyFace);
  MPI_Type_commit(&xzFace);
  MPI_Type_commit(&yzFace);

  if(rank == 0){
    
    /* Variable declaration/allocation. */
    full_cube = (double*)calloc(nx * ny * nz, sizeof(double));
    
    /* Initialize array(s). */
    // init_array(nx, full_cube);
    my_init_array(nx, ny, nz, full_cube);

    /* Start timer. */
    polybench_start_instruments;
  } 

  int *sendcounts = malloc(sizeof(int)*p);
  int *displs = malloc(sizeof(int)*p);
  int tmp_coords[3];
  for(int i = 0; i < p; i++){
    sendcounts[i] = 1;

    MPI_Cart_coords(topocomm, i, 3, tmp_coords);
    displs[i] = ((tmp_coords[0]*sx) * (ny) * (nz)) + ((tmp_coords[1]*sy) * (nz)) + (tmp_coords[2]*sz);
  }
  

  MPI_Scatterv(full_cube, sendcounts, displs, subCubeSend, IDX(A, 1, 1, 1, ssy, ssz), 1, subCubeRcv, 0, MPI_COMM_WORLD );


  for (int t = 1; t <= tsteps*2; t++) {


    const int num_requsts = 2 * 6;
    MPI_Request requests[num_requsts];
    MPI_Isend(IDX(A, 1    , 1, 1, ssy, ssz), 1, yzFace, front,  0, MPI_COMM_WORLD, requests + 0);
    MPI_Isend(IDX(A, ssx-2, 1, 1, ssy, ssz), 1, yzFace, back ,  0, MPI_COMM_WORLD, requests + 1);
    
    MPI_Isend(IDX(A, 1, 1, 1    , ssy, ssz), 1, xyFace, top,    0, MPI_COMM_WORLD, requests + 2);
    MPI_Isend(IDX(A, 1, 1, ssz-2, ssy, ssz), 1, xyFace, bottom, 0, MPI_COMM_WORLD, requests + 3);
    
    MPI_Isend(IDX(A, 1, 1    , 1, ssy, ssz), 1, xzFace, left,   0, MPI_COMM_WORLD, requests + 4);
    MPI_Isend(IDX(A, 1, ssy-2, 1, ssy, ssz), 1, xzFace, right,  0, MPI_COMM_WORLD, requests + 5);



    MPI_Irecv(IDX(A, 0    , 1, 1, ssy, ssz), 1, yzFace, front,  0, MPI_COMM_WORLD, requests + 6);
    MPI_Irecv(IDX(A, ssx-1, 1, 1, ssy, ssz), 1, yzFace, back,   0, MPI_COMM_WORLD, requests + 7);
    
    MPI_Irecv(IDX(A, 1, 1, 0    , ssy, ssz), 1, xyFace, top,    0, MPI_COMM_WORLD, requests + 8);
    MPI_Irecv(IDX(A, 1, 1, ssz-1, ssy, ssz), 1, xyFace, bottom, 0, MPI_COMM_WORLD, requests + 9);
    
    MPI_Irecv(IDX(A, 1, 0, 1    , ssy, ssz), 1, xzFace, left,   0, MPI_COMM_WORLD, requests + 10);
    MPI_Irecv(IDX(A, 1, ssy-1, 1, ssy, ssz), 1, xzFace, right,  0, MPI_COMM_WORLD, requests + 11);

    MPI_Waitall(num_requsts, requests, MPI_STATUS_IGNORE);


    compute_step_kernel_heat_3d(ssx, ssy, ssz, A, B);
  
  
    double *tmp = A;
    A = B;
    B = tmp;


  }


 /*
  for(int i = 0; i< p; i++){
 
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == i){
      fprintf(stderr, "\n\n==============\nScattered %d:\n==============\n", i);
      my_print_array(ssx, ssy, ssz, A);
      fprintf(stderr, "==============\n     %4d\n%4d %4d %4d\n     %4d\n==============\n", top, left, rank, right, bottom);
      fprintf(stderr, "==============\n     %4d\n%4d %4d %4d\n     %4d\n==============\n", top, front, rank, back, bottom);
    }
  }*/





  for(int i = 0; i < p; i++){
    sendcounts[i] = 1;

    MPI_Cart_coords(topocomm, i, 3, tmp_coords);
    displs[i] = ((tmp_coords[0]*sx) * (ny) * (nz)) + ((tmp_coords[1]*sy) * (nz)) + (tmp_coords[2]*sz);
  }


  MPI_Gatherv(IDX(A, 1, 1, 1, ssy, ssz), 1, subCubeRcv, full_cube, sendcounts, displs, subCubeSend, 0, MPI_COMM_WORLD);




  /* Be clean. */
  if(rank == 0){
    polybench_stop_instruments;
    polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
    polybench_prevent_dce(print_array(nx, full_cube));

    free(full_cube);
  }

  free(sendcounts);
  free(displs);

  
  free(A);
  free(B);
  MPI_Type_free(&subCubeSend);
  MPI_Type_free(&subCubeRcv);

  MPI_Type_free(&xyFace);
  MPI_Type_free(&xzFace);
  MPI_Type_free(&yzFace);

  MPI_Finalize();



  return 0;
}
