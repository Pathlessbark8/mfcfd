module parameter_mod
!
!
	implicit none
!
!
!	integer, parameter :: max_points=153600
!
	integer, parameter :: max_iters=500
!
!	real*8, parameter :: Mach=1.2d0
!	real*8, parameter :: aoa = 7.0d0
	real*8, parameter :: Mach=0.63d0
	real*8, parameter :: aoa = 2.0d0
!
!
	real*8, parameter :: rho_inf = 1.0d0
	real*8, parameter :: pr_inf = 1.0d0/1.4d0
	real*8, parameter :: gamma = 1.4d0	
!
!
!	The paramter theta 
!
	real*8, parameter :: pi=4.d0*datan(1.d0)				
	real*8, parameter :: theta = aoa*pi/180.d0
!
!
!
!	The parameter CFL is the CFL number for stability ..
!
	real*8, parameter :: CFL = 0.1d0
!
!	The parameter power is used to specify the weights 
!	in the LS formula for the derivatives ..
!	power = 0.0d0, -2.0d0, -4.0d0, -6.0d0 ..
!	For example, power = -2.0 implies that
!	power = -2.0 => weights = 1/d^2
!	power = -4.0 => weights = 1/d^4
!
	real*8, parameter :: power = 0.0d0 
!
!
!	limiter_flag = 1 => venkatakrishnan limiter
!	limiter_flag = 2 => min-max limiter 	
!
	integer, parameter :: limiter_flag = 1
	real*8, parameter :: VL_CONST = 150.d0 ! Venkatakrishnan limiter constant ..
!
!
!
!	flag for initial conditions. 
!	if flag = 0 => free stream conditions ..
!	if flag = 1 => read initial conditions from a file ..
!
	integer, parameter :: initial_conditions_flag = 0
!	
!	Interior points normal flag ..
!	If flag is zero => nx = 0.0 and ny = 1.0
!
	integer, parameter :: interior_points_normal_flag = 0
!	
!
!	Starting and end indices of various types of points ..
!	
!
!	integer, parameter :: interior_points = 152321 
!
!	integer, parameter :: outer_points = 640
!
!	integer, parameter :: wall_points = 639
!
	integer, parameter :: shapes = 1
!	integer, parameter :: shape_points(1) = 640	
	integer, parameter :: max_shape_points = 1280
!
!	Note: max_shape_points is the maximum value of all 
!	shape points (Maximum of shape_points_01, shape_points_02, ...)
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
	
