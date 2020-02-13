module parameter_mod

    implicit none

    real*8 :: Mach
    real*8 :: aoa
    real*8 :: theta
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
    real*8, parameter :: gamma = 1.4d0
	real*8, parameter :: pi=4.d0*datan(1.d0)
    real*8, dimension(4) :: q_init, q_inf

contains

    subroutine setup_case_parameters()
        implicit none

            !calculate theta
            theta = aoa*pi/180.d0

            !Setup initial conditions
            q_init(1) = rho_inf
            q_init(2) = Mach*dcos(theta)
            q_init(3) = Mach*dsin(theta)
            q_init(4) = pr_inf

            !Setup free stream conditions
            q_inf(1) = rho_inf
            q_inf(2) = Mach*dcos(theta)
            q_inf(3) = Mach*dsin(theta)
            q_inf(4) = pr_inf

    end subroutine

end module parameter_mod
