# Installation
The Meshfree Solver uses Python 3.x and cmake to setup the environment and compile the executable. 

## Prerequisites
Meshfree Solver requires

### Core Dependencies
* Python 3.x (for self-installer)
* HDF5 (For reading **.h5** files)

Apart from these dependencies the solver might require a number of transitive dependencies depending on the solver being compiled.

### Serial Solver Dependency
* PETSc, MPI (OpenMPI or MPICH)

### Tangent Solver Dependency
* PETSc, MPI (OpenMPI or MPICH)

### Adjoint Solver Dependency
* PETSc, MPI (OpenMPI or MPICH)

### CUDA Solver Dependency
* CUDA
* NVIDIA HPC SDK

## Package Outline 

There are four versions of the solver present in the four main directories i.e
* `src_mpi_serial` contains the code for the MPI version of the primal solver. 
* `src_mpi_cuda` contains the code for the CUDA version of the primal solver.
* `src_mpi_forward` contains the code for the MPI version of the differentiated tangent solver.
* `src_mpi_adjoint` contains the code for the MPI version of the differentiated adjoint solver.

In addition to these, we have
* `src_mpi_clean` contains the primal solver code (basic) from which the differentiated tangent and adjoint solver codes are created. 
* `install` contains the Python script to generate the executable for main verseions of the solvers.
* `examples` contains the subdirectories for the versions of the code containing all the required files and inputs via which the executable will run.

## Building
Run the following command from the `install` directory:

```
./install.py [-h] --mfcfd {cuda,serial,tangent,adjoint}
                  [--extra EXTRA_FLAGS]
```

Please specify the either cuda, serial, tangent or adjoint for the flag `mfcfd`.  
For example: `python3 install.py --mfcfd cuda`

After building, an executable `execname` will be generated.

## Environment Setup

Move the generated executable `execname`, created in the install directory to the required solver directory in `mfcfd/examples`   
For example, if the solver for adjoint has been compiled, then move the generated `execname` to inside the `mfcfd/examples/mpi_adjoint`

Move to the targeted subdirectory in `mfcfd/examples` and
* If the subdirectory is `cuda`, rename `input.mnl.example` to `input.nml`
* If the subdirectory is `mpi_adjoint`, `mpi_forward` or `mpi_serial`, rename `case.in.example` to `case.in` 

Finally in the same directory that the excutable has been shifted to, create a subdirectory called `point` and store the required HDF5 grid file to be run in that directory.  
For example if we want to do this for adjoint solver, `mkdir point` such that the path is `mfcfd/examples/mpi_forward/point` and store `point.h5` in this newly created subdirectory such that its path will be `mfcfd/examples/mpi_forward/point/point.h5`


## Execution

* For `serial`, `adjoint` and `tangent` versions of the solver, the executable is executed via the command
```
mpirun -np x ./execname
# X is the number of CPU cores on which the  solver will run
```

* For `cuda` version of the solver, the executable is executed via the command
```
./execname
```

## Output Files

*Here xxxx is the partitioned output given by a specific CPU core.  
For example: 0023 means the ouput was from the 24th core*

`serial`, `tangent` and `adjoint` solver outputs can be found in
* `cp` folder as `cp-filexxxx` which contains the cp output of the wall points in the form
```
point%flag_2(m) | point%x(m) | cp
# m is a wall point index. 
# If the file is blank then it did not have any wall points
```
* `solution` folder as `solxxxx.dat` which contains the solution output of all points in the form
```
point%original_id(i) | point%flag_1(i) | point%flag_2(i) | point%x(i)
| point%y(i) | point%prim(1,i) | point%prim(2,i) | point%prim(3,i) | point%prim(4,i)
```
* The adjoint solver in addition to the above files also generates the adjoint differentiated values of all points in the `sensitivity` folder as `sensitivity-xxxx.dat` in the form
```
point%original_id(i) | pointb%x(i) | pointb%y(i)
``` 

`cuda` solver outputs can be found in
* `cp-file` contains the cp output of the wall points in the form
```
point%flag_2(m) | point%x(m) | cp
```
* `output.dat` ontains the solution output of all points in the form
```
point%x(i) | point%y(i) | point%flag_1(i) | point%qtdepth(i) | point%prim(1,i) | point%prim(2,i) | point%prim(3,i) | point%prim(4,i) | mach_number | point%entropy(i) | point%sensor(i)
```

## Additional Information

* Major parameters like the `number of iterations`, `cfl`, `venkat_limter constant` can be changes in the required `input.nml` or `case.in` files 
* Most of the data structures and global variables are declared in the `data_structure.F90` or similarly named files
* Most of the file reading and associated memory alllocation subroutines can be found in the `point_preprocessor.F90` files.
* Most of the output related subroutines can be found in `post_processing.F90` files.
* Primary file for the solvers is the `meshfree_solver.F90` file.
* The actual computational subroutines can be found or are called in `q_lskum.F90` or similar named files and `fpi_solver.F90` files. 
* MPI/PETSc related communication calls and structures can be found in `petsc_data_structure.F90`
