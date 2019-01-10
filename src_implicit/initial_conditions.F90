module initial_conditions_mod

        use data_structure_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i

                if(initial_conditions_flag .eq. 0) then

                                do k=1, max_points
                                        point%prim(1,k) = rho_inf
                                        point%prim(2,k) = Mach*dcos(theta)
                                        point%prim(3,k) = Mach*dsin(theta)
                                        point%prim(4,k) = pr_inf
                                enddo

                endif

        end subroutine

end module initial_conditions_mod
