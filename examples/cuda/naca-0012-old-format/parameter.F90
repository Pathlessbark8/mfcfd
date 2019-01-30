module parameter_mod
#include "limiter.h"
        use cudafor

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
        real*8 :: power
!
        real*8 :: VL_CONST  ! Venkatakrishnan limiter constant ..

!       Interior points normal flag ..
!       If flag is zero => nx = 0.0 and ny = 1.0
!
        integer,parameter :: interior_points_normal_flag=0

!       Restart solution parameter
        integer :: solution_restart

!       old format tag
        integer,parameter :: old_format=1

!       First order flag
!        real*8,parameter, constant :: fo_flag=1.0

!       No of shapes
        integer,parameter :: shapes = 1

!       Block input
        integer :: blockx, blocky, blockz

        namelist / input_parameters /   &
                             cfl, &
                       max_iters, &
                          blockx, &
                          blocky, &
                          blockz, &
                          vl_const, &
                          power, &
                          solution_restart
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

        subroutine readnml()
      
                implicit none
                integer :: os

        
                open(unit=10,file='input.nml',form='formatted',status='old',iostat=os)
                read(unit=10,nml=input_parameters)
      
                close(10)
                write(*,*) '%%%%%%%%%%%%%-Input parameters-%%%%%%%%%%%%'
                write(*,*) 'Mach:', mach, 'aoa:', aoa
#ifdef VENKAT
                write(*,*) 'limiter:', 'venkat'
#endif
#ifdef MINMAX
                write(*,*) 'limiter:', 'minmax'
#endif
                write(*,*)
                write(*,*) '%%%%%%%%%%%%%%-Nml file info-%%%%%%%%%%%%%%'
                write(*,nml=input_parameters)
                write(*,*)
                



end subroutine


end module parameter_mod
