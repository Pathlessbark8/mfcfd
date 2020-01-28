program meshfree_solver

        use parameter_mod
        use data_structure_mod_diff
        use point_preprocessor_mod
        use q_lskum_mod_diff
        use q_lskum_mod_chkpts_diff
        ! use compute_force_coeffs_mod
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
        call allocate_soln_b()
!       Initialize Petsc vectors
        ! write(*,*) 'Number of points:         ', max_points
        write(*,*)

!       Assign the initial conditions for the primitive variables ..

        call initial_conditions()
        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
        write(*,*)

!	Primal fixed point iterative solver ..

        ! runtime = MPI_Wtime()
        ! call q_lskum_b()
        if(ad_mode == 0) then ! Black Box approach
                ! if(rank == 0) then
                        write(*,*)'%%%%%%%%%-Using Black box approach-%%%%%%%%'
                        write(*,*)
                ! end if
                call q_lskum_b()
        else
                ! if(rank == 0) then
                        write(*,*)'%%%%%%%-Using Checkpointing approach-%%%%%%'
                        write(*,*)
                ! end if
                call q_lskum_chkpts_b()
        end if

        ! call q_lskum_chkpts_b()
        ! runtime = MPI_Wtime() - runtime

!       Save solution one last time
        write(*,*) '%%%%%%%%%%%-Primal-%%%%%%%%%%%'
        call print_primal_output()

!       destroy petsc vectors and deallocate point/solution vectors
        write(*,*) '%%%%%%%%%%%-Dealloc Soln-%%%%%%%%%%%'
        call deallocate_soln()
        write(*,*) '%%%%%%%%%%%-Dealloc Backwards Soln-%%%%%%%%%%%'
        call deallocate_soln_b()
        call dealloc_points()

        ! totaltime = MPI_Wtime() - totaltime

        write(*,*)
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',totaltime,'seconds'

!       stop petsc
        ! call MPI_Finalize(ierr)

end program meshfree_solver
