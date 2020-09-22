# Installation
The Meshfree Solver uses Python 3.x and cmake to setup the environment and compile the executable. 

## Prerequisites
Meshfree Solver requires

### Core Dependencies
* Python 3.x (for self-installer)
* HDF5 (For reading **.h5** files)

Apart from these dependencies the solver might require a number of transitive dependencies depending on the solver being compiled.

### Serial Solver Dependency
* PETSc

### Tangent Solver Dependency
* PETSc

### Adjoint Solver Dependency
* PETSc

### CUDA Solver Dependency
* CUDA
* NVIDIA HPC SDK

## Building
Run the following command from the `install` directory:

```
./install.py [-h] --mfcfd {cuda,serial,tangent,adjoint}
                  [--extra EXTRA_FLAGS]
```

Please specify the either cuda, serial, tangent or adjoint for the flag `mfcfd`.

