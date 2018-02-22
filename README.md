# FreeFLOW

FreeFLOW is a 3D, parallel Meshfree solver integrated with Adjoint based Algorithmic differentiation model. The code is made parallel using the data types availabe in the High Performance Computing tool [PETSc](https://www.mcs.anl.gov/petsc/).

#Prerequisites 

1. Download spack, [here](https://spack.io/)
2. Install PETSc.
3. Set the PETSc path to `$PETSC_DIR`.

#For Users

1. To clone the repository firstly install git on your system.
2. Setup an SSH key in your computer and further link it your bitbucket account, see [here](https://confluence.atlassian.com/bitbucket/set-up-an-ssh-key-728138079.html).
3. Later clone the repository using the following command in your terminal,
`git clone git@bitbucket.org:srikanthcs05/freeflow.git`.
4. Spend some time on bit bucket tutorials,[here](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud).

#Running the program

1. First, set the path to the cloned FreeFLOW directory to `$FF_DIR`.
2. Go to examples directory and choose a test case.
3. Make necessary changes in the paramter.f90 folder.
4. Run `make all`.
5. An executable `execname` will be generated. Run it.





#Developers

1. [Dr. N. Anil](http://universe.bits-pilani.ac.in/hyderabad/nanil/Profile)
2. Srikanth Chowdadenahalli Sathyanarayana
3. Shubham Ranjan

