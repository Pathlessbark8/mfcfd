module initial_conditions_mod

    use data_structure_mod_diff
    use petsc_data_structure_mod
    use parameter_mod

contains

    subroutine initial_conditions()

        implicit none
        
        integer :: k,i

        if(restart == 0) then

            call setup_case_parameters()

            do k=1, max_points
                point%prim(1,k) = q_init(1)
                point%prim(2,k) = q_init(2)
                point%prim(3,k) = q_init(3)
                point%prim(4,k) = q_init(4)
            enddo
        
        elseif(restart == 1) then
    
            if(rank == 0) then
                write(*,*)'%%%%%%%%%%%%-Using restart file-%%%%%%%%%%%'
                write(*,*)
            end if

            call setup_case_parameters()
            
            call restart_sol()

            call update_begin_prim_ghost()
            call update_end_prim_ghost()
        
        endif

    end subroutine

    subroutine restart_sol()
 
        implicit none

        integer :: i, dummy
        character(len=64) :: sfile
        character(len=10) :: itos

        if(proc==1) then
            sfile = 'restart/sol.dat'
        else
            sfile = 'restart/'//'sol-'//trim(itos(4,rank))//'.dat'
        end if

        OPEN(UNIT=515,FILE=trim(sfile),form='formatted', action="read")

        read(515,*)dummy, itr, res_old

        do i = 1, local_points
            read(515,*)dummy, dummy, dummy, dummy,&
                dummy, point%prim(1,i), point%prim(2,i), point%prim(3,i),&
                point%prim(4,i)
        end do

        close(515)

    end subroutine

end module initial_conditions_mod
