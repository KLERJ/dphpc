# how to use benchmark.py

## general
* the script looks in the `src/` folder for executables
* make sure the executable you pass contains either "deriche", "seidel-2d" or "heat-3d"
* you need to generate the reference solutions using the `generate_refs.sh` scripts beforehand

## flags
* `-s` or `--dataset-size` is one of MINI, SMALL, MEDIUM, LARGE or EXTRALARGE
* `-n` or `--runs` is the number of runs to average over
* `-np` or `--procs` number of processors to use
* `-t` or `--test` whether to test the implementation. If absent, the implementation is benchmarked
* `-d` or `--diff` whether to show the diff of the output files. If absent, you just get binary feeback
if the outputs match

## testing the correctness of an implementation
* run it like `./benchmark.py seidel-2d_omp --dataset-size MINI --test`
* this will tell you whether the output is the same as the one from the reference
solution
* if you would like to see the diff, add the `--diff` or `-d` flag to the above command

## benchmarking the implementation
* run it like `./benchmark.py seidel-2d_omp --dataset-size MINI --runs 5`
* this will average over 5 runs and write the runtimes to a `.csv` file
