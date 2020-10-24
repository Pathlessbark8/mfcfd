module stagnation_values_mod

    use device_data_structure_mod

    contains

    attributes(global) subroutine objective_function_J(prim_d, cost_func)

        implicit none

        ! device variables
        real*8 :: prim_d(:,:)
        real*8 :: cost_func(:)

        ! local variables
        real*8 :: p0_inf, gammaPower, p0, p0_sum, constant, angle, mach_t
        real*8 :: prim(4)
        integer :: i


        i = (blockIdx%x-1)* blockDim%x + threadIdx%x

        if(i > mp_d) return

        gammaPower = gamma_d/(gamma_d-1)
        p0_inf = pr_inf*((1 + ((gamma_d - 1)/2)*mach_d*mach_d) ** gammaPower)

        constant = 1/(p0_inf**2 * mp_d)
        p0_sum = 0.d0
        prim = prim_d(:,i)
        angle = sqrt(gamma_d * prim(4)/ prim(1))
        mach_t = sqrt(prim(2)**2 + prim(3)**2)/angle
        p0 = prim(4)*((1 + ((gamma_d - 1)/2)*mach_t*mach_t) ** gammaPower)
        cost_func(i) = ((p0_inf - p0) ** 2) * constant

    end subroutine
    
end module stagnation_values_mod
