# An Implicitly Parallel Meshfree Solver in Fortran

## Dependencies
1. [PETSc](https://www.mcs.anl.gov/petsc/) 
2. GNU make
3. MPI (OpenMPI 3.x / MPICH 3.x)
4. GCC 7.x and above


## Optional Dependencies
1. Python 3.x


# Installation

## For PAW-ATM 2020
All results and benchmarks were computed using the `mpi-second-order` branch.

### 40M Grid File
We have provided a 40 Million fine grid `partGrid40M_unpartitioned` which can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1iqPZOxj0UBDS3u6mv-CAHdOgbhXKAIYN).

Reviewers are requested to use the [mfpre](https://github.com/Nischay-Pro/mfpre) partitioner with `--quadtree` flag along with their specified number of partitions they want.

### Example Script
An example make file can be accessed from `examples/mpi_serial`. Reviewers can use this folder to compile and execute the solver.

### Compilation
1. Copy `src_mpi_serial/parameter.F90` file to the current folder containing the make file.
2. Set the environment variable MF_DIR. 
 ```
 export MF_DIR=<path to the top-level directory where the respository is downloaded>
 # Example
 # export MF_DIR=/home/nischay/Git/mfcfd
 ```
3. Load PETSc, GCC and MPI modules for compilation.
4. Run `make`
5. In the current folder a `execname` executable will be created. 

### Configuration of the test case.
1. Run the partitioner on the grid file provided in Google Drive.
2. A folder called `point` will be generated in the directory of the partitioner.
3. Using `mv` shift `point` to the folder containing the compiled executable.
4. Edit the provided `case.in` file with the following configuration.

```
-cfl 0.2
-max_iters 1000
-tscheme ssprk43
-format_file quadtree
-vl_const 50
-inner_iterations 3
-mach 0.85
-aoa 1
```

### Execution
1. Launch the `execname` executable with the number of partitions done on the grid. 
2. For example if you have 16 partitions use `mpirun -np 16 ./execname`.