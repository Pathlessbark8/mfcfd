program meshfree_solver

#include <petsc/finclude/petscsys.h>


        use petscsys
	use parameter_mod
	use data_structure_mod
	use point_preprocessor_mod
        use initial_conditions_mod
	use q_lskum_mod
!	use post_processing_mod


	implicit none

        PetscErrorCode  :: ierr


!
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        if(ierr /= 0) stop "Unable to initialize PETSc"
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, proc, ierr); CHKERRQ(ierr)
!
!	Reading the input data ..
!
	call read_input_point_data()
!
!	Assign the initial conditions for the primitive variables ..	
!
	call initial_conditions()

         
        if (proc > 1) call init_petsc()
        print*,'working till here'
!
!	Primal fixed point iterative solver ..
!
	call q_lskum()
!
!	Printing the output (post-processing) ..
!
!	call print_primal_output()

        if (proc>1) call dest_petsc()

        call PetscFinalize(ierr); CHKERRQ(ierr)
	
end program meshfree_solver
