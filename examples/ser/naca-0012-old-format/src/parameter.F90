module parameter_mod

        implicit none

        real*8, parameter :: Mach = 0.63d0
        real*8, parameter :: aoa = 2.0d0
        real*8 :: theta
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
        real*8, parameter :: gamma = 1.4d0
	real*8, parameter :: pi=4.d0*datan(1.d0)
        real*8, dimension(4) :: q_init, q_inf

        !The parameter CFL is the CFL number for stability ..
        real*8 :: CFL

        integer :: max_iters
!
!       The parameter power is used to specify the weights 
!       in the LS formula for the derivatives ..
!       power = 0.0d0, -2.0d0, -4.0d0, -6.0d0 ..
!       For example, power = -2.0 implies that
!       power = -2.0 => weights = 1/d^2
!       power = -4.0 => weights = 1/d^4
!
        real*8,parameter :: power=0.0d0
!
!
!       limiter_flag = 1 => venkatakrishnan limiter
!       limiter_flag = 2 => min-max limiter     
!
        integer,parameter :: limiter_flag=1
        real*8,parameter :: VL_CONST=150.d0  ! Venkatakrishnan limiter constant ..

        integer,parameter :: initial_conditions_flag=0

!       Interior points normal flag ..
!       If flag is zero => nx = 0.0 and ny = 1.0
!
        integer,parameter :: interior_points_normal_flag=0

!       Restart solution parameter
        integer :: solution_restart

!       solution save parameter
        integer :: nsave

!       Time scheme, explicit or implicit
        integer :: tscheme

!       old format tag
        integer,parameter :: old_format=1

!       First order flag
        real*8,parameter :: fo_flag=1.0

!       Objective function
        real*8 :: Cl_flag, Cd_flag, Cm_flag, Cl_Cd_flag, ent_flag, ens_flag
        integer,parameter :: obj_flag=0

!       No of shapes
        integer,parameter :: shapes=1


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

        subroutine readinp()
      
                implicit none
      
                open(unit=11,file='input.dat',form='formatted')
                read(11,*)
                read(11,*) max_iters, cfl
      
                close(11)

end subroutine


end module parameter_mod
