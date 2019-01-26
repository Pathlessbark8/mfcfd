module initial_conditions_mod

        use data_structure_mod
        use parameter_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i, dummy
                real*8 :: dummyx, dummyy


                call setup_case_parameters()

                if(solution_restart .eq. 0) then

                        do k=1, max_points
                                point%prim(1,k) = q_init(1)
                                point%prim(2,k) = q_init(2)
                                point%prim(3,k) = q_init(3)
                                point%prim(4,k) = q_init(4)
                        enddo
                else

                        write(*,*) '%%%%%%%%%%%%-Using Restart file-%%%%%%%%%%%%'

                        open(unit=105, file="restart.dat", form='formatted', action="read")

                        do k=1, max_points
                               read(105,*) dummy, dummy, dummy, dummyx, dummyy, point%prim(1,k), point%prim(2,k), point%prim(3,k), point%prim(4,k)
                        end do

                        close(unit=105)

                endif

        end subroutine

end module initial_conditions_mod