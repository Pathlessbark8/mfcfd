program meshfree_solver

        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use q_lskum_mod
        use compute_force_coeffs_mod

        implicit none
        real*8  :: totaltime,runtime
        integer :: ierr

        call MPI_INIT(ierr)
        if(ierr /= 0) stop "Unable to initialize MPI"
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, proc, ierr)
        
        if(rank==0) then
                call execute_command_line('mkdir -p solution')
                call execute_command_line('mkdir -p cp')
        end if
        
        totaltime = MPI_Wtime()

        if(rank == 0) then
                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-MPI Meshfree Code-%%%%%%%%%%%'
                write(*,*)
        end if

!       Read the case file

        if(rank == 0) then
                write(*,*)'%%%%%%%%%%%-Reading the case file-%%%%%%%%%'
                write(*,*)
        end if
        
        call readcase()

!	Reading the input data ..

        if(rank == 0) then
                write(*,*)'%%%%%%%%%%%%-Reading point data-%%%%%%%%%%%'
                write(*,*)
        end if
        
        call read_input_point_data()
        
!       Allocate solution variables

        call allocate_soln()

!       Initialize Petsc vectors

        if(proc == 1) plen = max_points
        if(rank == 0) then
                write(*,*) 'Number of points:         ', plen
                write(*,*)
        end if

!	Primal fixed point iterative solver ..
        
        runtime = MPI_Wtime()
        call q_lskum()
        runtime = MPI_Wtime() - runtime

!       Save solution one last time
        call print_primal_output()

!       destroy petsc vectors and deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()

        totaltime = MPI_Wtime() - totaltime

        if(rank == 0) then
                write(*,*)
                write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
                write(*,*) 'Run time:  ',runtime,'seconds'
                write(*,*) 'Total time:',totaltime,'seconds'
        end if

!       stop petsc
        call MPI_Finalize(ierr)

end program meshfree_solver
