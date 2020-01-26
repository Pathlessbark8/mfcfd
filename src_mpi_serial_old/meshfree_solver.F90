program meshfree_solver
        
        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use q_lskum_mod
        ! use compute_force_coeffs_mod
        use file_ops_mod
        use initial_conditions_mod
        use post_processing_mod

        implicit none
        INCLUDE 'mpif.h'

        real*8  :: totaltime,runtime
        integer :: ierr

        call MPI_INIT(ierr)

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, proc, ierr)
        if(rank==0) then
                call execute_command_line('mkdir -p solution')
                call execute_command_line('mkdir -p cp')
        end if     
        ! totaltime = MPI_Wtime()

        
        if(rank == 0) then
                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-MPI Meshfree Code-%%%%%%%%%%%'
                write(*,*)
        end if

!	Reading the input data ..

        if(rank == 0) then
                write(*,*)'%%%%%%%%%%%-Reading the case file-%%%%%%%%%'
                write(*,*)
        end if

        call readnml()

        write(*,*)'%%%%%%%%%%%%-Reading point data-%%%%%%%%%%%'
        write(*,*)

        call read_input_point_data()

!       Allocate solution variables

        call allocate_soln()
        ! call allocate_soln_d()
!       Initialize Petsc vectors
        write(*,*) 'Number of points:         ', max_points
        write(*,*)

!       Assign the initial conditions for the primitive variables ..

        call initial_conditions()
        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
        write(*,*)

!	Primal fixed point iterative solver ..

        ! runtime = MPI_Wtime()
        ! call q_lskum()
        ! runtime = MPI_Wtime() - runtime

!       Save solution one last time
        call print_primal_output()

!       destroy petsc vectors and deallocate point/solution vectors
        call deallocate_soln()
        ! call deallocate_soln_d()
        call dealloc_points()

        ! totaltime = MPI_Wtime() - totaltime

        write(*,*)
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',totaltime,'seconds'

!       stop petsc
        call MPI_Finalize(ierr)

end program meshfree_solver
