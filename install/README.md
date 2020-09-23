# Installation
The Meshfree Solver uses Python 3.x and cmake to setup the environment and compile the executable. 

## Prerequisites
Meshfree Solver requires

### Core Dependencies
* Python 3.x (for self-installer)
* HDF5 (For reading **.h5** files)

Apart from these dependencies the solver might require a number of transitive dependencies depending on the solver being compiled.

### Serial Solver Dependency
* PETSc, MPI

### Tangent Solver Dependency
* PETSc, MPI

### Adjoint Solver Dependency
* PETSc, MPI

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
