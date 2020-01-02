program meshfree_solver

        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use q_lskum_mod
        use compute_force_coeffs_mod
        use file_ops_mod
        use initial_conditions_mod
        use post_processing_mod

        implicit none
        real*8  :: totaltime,runtime
        integer :: ierr
        
        call execute_command_line('mkdir -p solution')
        call execute_command_line('mkdir -p cp')
        
        ! totaltime = MPI_Wtime()

        write(*,*)
        write(*,*)'%%%%%%%%%%%%%-Serial Meshfree Code-%%%%%%%%%%%'
        write(*,*)

!	Reading the input data ..

        call readnml()

        write(*,*)'%%%%%%%%%%%%-Reading point data-%%%%%%%%%%%'
        write(*,*)
        
        call read_input_point_data()
        
!       Allocate solution variables

        call allocate_soln()

!       Initialize Petsc vectors
        write(*,*) 'Number of points:         ', max_points
        write(*,*)

!       Assign the initial conditions for the primitive variables ..	

        call initial_conditions()
        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
        write(*,*)

!	Primal fixed point iterative solver ..
        
        ! runtime = MPI_Wtime()
        call q_lskum()
        ! runtime = MPI_Wtime() - runtime

!       Save solution one last time
        call print_primal_output()

!       destroy petsc vectors and deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()

        ! totaltime = MPI_Wtime() - totaltime

        write(*,*)
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',totaltime,'seconds'

!       stop petsc
        ! call MPI_Finalize(ierr)

end program meshfree_solver
