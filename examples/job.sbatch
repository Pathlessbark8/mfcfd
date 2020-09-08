#!/bin/bash

#SBATCH -p gpu          # Submit to 'stdmemq' Partitiion or queue
#SBATCH -J FortranGPU          # Name the job as 'MPItest'
#SBATCH -o output/GPUtest-%j.out   # Write the standard output to file named 'jMPItest-<job_number>.out'
#SBATCH -e output/GPUtest-%j.err   # Write the standard error to file named 'jMPItest-<job_number>.err'
#SBATCH -t 0-12:00:00        # Run for a maximum time of 0 days, 12 hours, 00 mins, 00 secs
#SBATCH --ntasks-per-node=1 # Request n cores or task per node
#SBATCH --mem-per-cpu=4000   # Request 4000MB (4GB) RAM per core
#SBATCH --exclusive
#SBATCH --gres=gpu:1
#SBATCH --reservation=meshfree_217

TMPDIR=/home/meshfree/temp_files/
module load nvhpc/20.7
#srun --gres=gpu:1 ./execname
# srun --gres=gpu:1 /home/opt_ohpc_pub/nvhpc-sdk/Linux_x86_64/20.7/compilers/bin/nv-nsight-cu-cli -f -o profile/run1_initial ./execname
srun --gres=gpu:1 nv-nsight-cu-cli -f -o profile/run2_initial_64threads_85reg ./execname
#srun --gres=gpu:1 nvprof ./execname
