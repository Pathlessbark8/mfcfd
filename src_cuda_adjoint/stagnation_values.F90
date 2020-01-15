module stagnation_values_mod

    use DATA_STRUCTURE_MOD_DIFF

    contains

        subroutine stagnation_pressure()

            implicit none

            integer :: i, indexMin, indexMax
            real*8 :: p0_inf, gammaPower, pMin, pMax, p0, mach_t, angle
            real*8 :: prim(4)

            gammaPower = gamma/(gamma-1)
            p0_inf = pr_inf*((1 + ((gamma - 1)/2)*mach*mach) ** gammaPower)

            do i=1, max_points
                prim = point%prim(:,i)
                angle = sqrt(gamma * prim(4)/ prim(1))
                mach_t = sqrt(prim(2)**2 + prim(3)**2)/angle
                p0 = prim(4)*((1 + ((gamma - 1)/2)*mach_t*mach_t) ** gammaPower)
                if(i == 1) then
                    pMin = p0
                    indexMin = i
                    indexMax = i
                    pMax = p0
                elseif(p0 < pMin) then
                    pMin = p0
                    indexMin = i
                elseif(p0 > pMax) then
                    pMax = p0
                    indexMax = i
                endif
            enddo

            write(*,*) "Stagnation values are ", pMin, " ", pMax, " ", pMin/p0_inf, " ", pMax/p0_inf," ", indexMin, " ", indexMax

        end subroutine

        subroutine objective_function_J()

            implicit none

            integer :: i
            real*8 :: p0_inf, gammaPower, p0, p0_sum, constant, angle, mach_t
            real*8 :: prim(4)

            gammaPower = gamma/(gamma-1)
            p0_inf = pr_inf*((1 + ((gamma - 1)/2)*mach*mach) ** gammaPower)

            constant = 1/(p0_inf**2 * max_points)
            p0_sum = 0

            do i=1, max_points
                prim = point%prim(:,i)
                angle = sqrt(gamma * prim(4)/ prim(1))
                mach_t = sqrt(prim(2)**2 + prim(3)**2)/angle
                p0 = prim(4)*((1 + ((gamma - 1)/2)*mach_t*mach_t) ** gammaPower)
                p0_sum = p0_sum + (p0_inf - p0) ** 2
            enddo

            cost_func = p0_sum * constant

            write(*,*) "Objective Function (J)", cost_func

        end subroutine


end module stagnation_values_mod
