module data_structure_mod
        
        use parameter_mod

        implicit none

        integer :: max_points=9600, local_points=9600
        integer :: wall_points=159,interior_points=9281,outer_points=160,shape_points=160

        type :: points


                integer, dimension(max_points) :: original_id
		real*8, dimension(max_points) :: x,y
                integer, dimension(max_points) :: left,right
                integer, dimension(max_points) :: flag_1 ! stores location of point
                integer, dimension(max_points) :: flag_2 ! stores shape point belongs to 
		real*8, dimension(max_points) :: nx,ny
                integer, dimension(max_points) :: nbhs
                integer, dimension(max_points,15) :: conn

		real*8, dimension(max_points) :: min_dist

                real*8, dimension(4,max_points) :: prim
                real*8, dimension(4,max_points) :: prim_old
		real*8, dimension(4,max_points) :: flux_res

                real*8, dimension(4,max_points) :: q
		real*8, dimension(2,4,max_points) :: dq
		real*8, dimension(2,4,max_points) :: qm

		real*8, dimension(max_points) :: entropy, vorticity, vorticity_sqr

                integer, dimension(max_points) :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer, dimension(max_points,15) :: xpos_conn, xneg_conn
                integer, dimension(max_points,15) :: ypos_conn, yneg_conn

                real*8, dimension(max_points)  :: delta
                
! Implicit data
                real*8, dimension(4,max_points) :: U_old
        end type points
 
        type(points) :: point

        save

        integer,dimension(wall_points) :: wall_points_index
        integer,dimension(wall_points) :: outer_points_index
        integer,dimension(wall_points) :: interior_points_index
        integer,dimension(wall_points) :: shape_points_index

        !iterations
        integer :: it


        real*8  :: res_old, res_new, residue, max_res
        real* 8 :: gsum_res_sqr,sum_res_sqr
        integer :: max_res_point
	real*8, dimension(shapes)  :: Cl, Cd, Cm, ClCd
	real*8  :: total_entropy, total_enstrophy
        integer :: plen

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
!       limiter_flag = 1 => venkatakrishnan limiter
!       limiter_flag = 2 => min-max limiter     
!
        integer :: limiter_flag
        real*8 :: VL_CONST  ! Venkatakrishnan limiter constant ..

        integer :: initial_conditions_flag

!       Interior points normal flag ..
!       If flag is zero => nx = 0.0 and ny = 1.0
!
        integer :: interior_points_normal_flag

!       Restart solution parameter
        integer :: solution_restart

!       solution save parameter
        integer :: nsave

!       First order flag
        real*8 :: fo_flag

!       Objective function
        real*8 :: Cl_flag, Cd_flag, Cm_flag, Cl_Cd_flag, ent_flag, ens_flag
        integer :: obj_flag

!       No of shapes
        integer :: shapes=1

end module data_structure_mod
