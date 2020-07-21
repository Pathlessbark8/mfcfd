module initial_conditions_mod

        use data_structure_mod
        use petsc_data_structure_mod
        use parameter_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i

                CALL UPDATE_BEGIN_PHI1_GHOST()
                CALL UPDATE_END_PHI1_GHOST()
                CALL UPDATE_BEGIN_PHI2_GHOST()
                CALL UPDATE_END_PHI2_GHOST()

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

                call voronoi_area_read()

        end subroutine

        subroutine restart_sol()
 
                implicit none

                integer :: i
                character(len=64) :: sfile
                character(len=10) :: itos
                real*8 :: dummy

                if(proc==1) then
                        sfile = 'restart/sol.dat'
                else
                        sfile = 'restart/'//'sol-'//trim(itos(4,rank))//'.dat'
                end if

                OPEN(UNIT=515,FILE=trim(sfile),form='formatted', action="read")

                read(515,*)dummy, itr, res_old

                do i = 1, local_points
                        read(515,*)dummy, dummy, dummy, dummy,&
                                point%prim(1,i), point%prim(2,i), point%prim(3,i),&
                                point%prim(4,i), dummy, dummy, dummy
                end do

                close(515)


        end subroutine

        subroutine voronoi_area_read()
 
            implicit none

            integer :: i
            character(len=64) :: sfile
            character(len=10) :: itos
            real*8 :: dummy
            
            if(proc==1) then
                    sfile = 'voronoi_values_split/voronoi.dat'
            else
                    sfile = 'voronoi_values_split/'//'voronoi'//trim(itos(4,rank))//'.dat'
            end if

            OPEN(UNIT=515,FILE=trim(sfile),form='formatted', action="read")

            do i = 1, local_points
                read(515,*) dummy, point%vor_area(i)
            end do

            close(515)

        end subroutine

end module initial_conditions_mod
