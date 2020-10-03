module q_lskum_mod

!	First written on 14.10.2016
!	updated on Dec 26, 2016
!	updated on Dec 29, 2016

    use data_structure_mod
    use point_normals_mod    
    use generate_connectivity_mod
    use fpi_solver_mod
    use initial_conditions_mod
    use ieee_arithmetic

contains

    subroutine q_lskum()

        implicit none

        integer :: i

        if(rank==0)OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
        ! if(point%original_id(1) == 1) then
        !     point%x(1) = point%x(1) + 1E-6
        ! end if
        ! Reminder to have x and /or y communicated if modified for ghost points
        call compute_normals()
        call generate_connectivity()

        if(rank == 0) then
            write(*,*)
            write(*,*)'%%%%-Normals and connectivity generated-%%%'
            write(*,*)
        end if
        ! Set U_old to U for first iteration
        do i=1,local_points
            point%U_old(1,i) = point%prim(1,i)
            point%U_old(2,i) = point%prim(1,i)*point%prim(2,i)
            point%U_old(3,i) = point%prim(1,i)*point%prim(3,i)
            point%U_old(4,i) = 2.5d0*point%prim(4,i) + 0.5d0*point%prim(1,i)*&
                &(point%prim(2,i)*point%prim(2,i) +&
                &point%prim(3,i)*point%prim(3,i))
        end do
        
        if(rank == 0) then
            write(*,*)'%%%%%%%%%%%%%-Iterations begin-%%%%%%%%%%%%'
            write(*,*)
        end if

        t = 0.0d0
        if(restart == 0)itr = 0
        
        do it = itr+1, itr+max_iters
            
            call fpi_solver(it)
            t = t + dtg
            if (rank==0) then
                if(timestep == 0) then
                    write(*,'(a12,i8,a15,e30.20)')'iterations:',it,'residue:',residue
                    write(301, *) it, residue
                elseif(timestep == 1) then
                    write(*,'(a12,i8,a15,e30.20)')'iterations:',it,'time:',t
                    write(301, *) it, t, dtg
                end if
                if(ieee_is_nan(residue))exit
            end if
        enddo
        
        CLOSE(UNIT=301)
        if(rank==0) then
            write(*,*) "Vector function is ", vector_cost_func
        end if
    end subroutine

end module q_lskum_mod
