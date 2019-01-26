module initial_conditions_mod

        use data_structure_mod
        use parameter_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i

                if(initial_conditions_flag .eq. 0) then

                        call setup_case_parameters()

                        do k=1, max_points
                                point%prim(1,k) = q_init(1)
                                point%prim(2,k) = q_init(2)
                                point%prim(3,k) = q_init(3)
                                point%prim(4,k) = q_init(4)
                        enddo

                endif

        end subroutine

end module initial_conditions_mod
