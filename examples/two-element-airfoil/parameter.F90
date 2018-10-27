module parameter_mod
!
!
	implicit none
!
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
	real*8, parameter :: gamma = 1.4d0	
!
!
	real*8, parameter :: pi=4.d0*datan(1.d0)				
!
!
!	Parameters for the adjoint code ..
!
!
!        integer, parameter :: obj_func_flag = 0

	real*8, parameter :: Cl_flag = 0.0d0
	real*8, parameter :: Cd_flag = 0.0d0
	real*8, parameter :: Cm_flag = 0.0d0
        real*8, parameter :: Cl_Cd_flag = 0.0d0
!
!	If Cl_flag = 1.0, then it is the objective function. 
! 	Otherwise if Cl_flag = 0.0, then it is not the 
!	objective function ..
!
!
!       Number of checkpoints ..
!
	integer, parameter :: ckpts = 350
!  
!
!  REAL :: DATAN
!  EXTERNAL DATAN
!	
!
end module parameter_mod
	
