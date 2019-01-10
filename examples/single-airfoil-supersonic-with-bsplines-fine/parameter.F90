module parameter_mod

        implicit none

! Set the freestream conditions

        real*8, parameter :: Mach = 1.2
        real*8, parameter :: aoa = 0.0d0
        real*8 :: theta
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
        real*8, parameter :: gamma = 1.4d0
	real*8, parameter :: pi=4.d0*datan(1.d0)
        real*8, dimension(4) :: q_init 

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

        end subroutine

end module parameter_mod
