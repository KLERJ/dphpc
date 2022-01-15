Design of Parallel and High-Performance Computing Project
======

## Fast Parallel Implementations of Deriche Edge Detector and Heat-3D Equation Solver using OpenMP and OpenMPI



### Compiling and running the code
```bash
# Build everything with default parameters
make

# Run benchmark with given parameters
mpirun -n 4 ./bin/heat-3d_mpi <N_SIZE> <T_STEPS> OUTPUT_DUMP BENCHMARK_DUMP_PATH

## The results are saved into benchmark_{0-3}.csv
# 0.085950,0.086273,0.069604,0.074362,0.086158,0.069257...
```
