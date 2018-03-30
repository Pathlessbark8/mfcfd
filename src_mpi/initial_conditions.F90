module initial_conditions_mod

        use data_structure_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i

                if(initial_conditions_flag .eq. 0) then

                                do k=1, max_points
                                        p%rho(k) = rho_inf
                                        p%u1(k) = Mach*dcos(theta)
                                        p%u2(k) = Mach*dsin(theta)
                                        p%pr(k) = pr_inf
                                enddo

                else if(initial_conditions_flag .eq. 1) then

                        !TODO Input file for partitions

                        OPEN(UNIT=100,FILE="./structured-grid-160-60/stored-solution",FORM="FORMATTED",STATUS="OLD",ACTION="READ")

                                do k=1, local_points
                                        read(100,*) i, p%rho(k), p%u1(k), p%u2(k), p%pr(k)
                                enddo

                                CLOSE(UNIT=100)
                endif

        end subroutine

end module initial_conditions_mod
