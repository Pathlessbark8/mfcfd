module parameter_mod

        implicit none

        real*8 :: Mach = 0.63d0
        real*8 :: aoa = 2.0d0
        real*8 :: aoad = 1.0d0
        real*8 :: theta, thetad
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
        real*8, parameter :: gamma = 1.4d0
	real*8, parameter :: pi=4.d0*datan(1.d0)
        real*8, dimension(4) :: q_init, q_inf
        real*8, dimension(4) :: q_initd, q_infd

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

        subroutine setup_case_parameters_d()
                
                implicit none
                !calculate theta
                thetad = pi*aoad/180.d0
                theta = aoa*pi/180.d0
                !Setup initial conditions
                q_init(1) = rho_inf
                q_initd = 0.0_8
                q_initd(2) = -(mach*thetad*DSIN(theta))
                q_init(2) = mach*DCOS(theta)
                q_initd(3) = mach*thetad*DCOS(theta)
                q_init(3) = mach*DSIN(theta)
                q_initd(4) = 0.0_8
                q_init(4) = pr_inf
                !Setup free stream conditions
                q_inf(1) = rho_inf
                q_infd = 0.0_8
                q_infd(2) = -(mach*thetad*DSIN(theta))
                q_inf(2) = mach*DCOS(theta)
                q_infd(3) = mach*thetad*DCOS(theta)
                q_inf(3) = mach*DSIN(theta)
                q_infd(4) = 0.0_8
                q_inf(4) = pr_inf
        end subroutine setup_case_parameters_d

end module parameter_mod
