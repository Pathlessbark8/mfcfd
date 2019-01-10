program meshfree_solver

#include <petsc/finclude/petscsys.h>


        use petscsys
        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use initial_conditions_mod
        use q_lskum_mod
        use compute_force_coeffs_mod


        implicit none
        real*8  :: totaltime,runtime
        PetscErrorCode  :: ierr


!
        call PetscInitialize('case.in', ierr)
        if(ierr /= 0) stop "Unable to initialize PETSc"
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, proc, ierr)
        if(rank==0) then
                call execute_command_line('mkdir -p solution')
                call execute_command_line('mkdir -p cp')
        end if

!
!	Reading the input data ..

        totaltime = MPI_Wtime()
!
        call read_input_point_data()

        !       Read the case file

        call readcase()

!       Allocate solution variables

        call allocate_soln()

!       Initialize Petsc vectors

        if(proc .ne. 1)call init_petsc()
        if(proc==1) plen = max_points
        if(rank==0) print*,'No of points:',plen

!	Assign the initial conditions for the primitive variables ..	

        call initial_conditions()
        if(rank==0)print*,'solution initialised'
       
!	Primal fixed point iterative solver ..
        
        runtime = MPI_Wtime()
        call q_lskum()
        runtime = MPI_Wtime() - runtime

!       Save solution one last time
        call print_primal_output()


!       destroy petsc vectors and deallocate point/solution vectors
        call dest_petsc()
        call deallocate_soln()
        call dealloc_points()

        totaltime = MPI_Wtime() - totaltime
        if(rank==0)print*,'Simulation completed!!!'
        if(rank==0)print*,'runtime:',runtime
        if(rank==0)print*,'totaltime:',totaltime
        
!       stop petsc
        call PetscFinalize(ierr)

end program meshfree_solver
