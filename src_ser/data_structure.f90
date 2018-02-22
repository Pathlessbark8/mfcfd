! flag_1 : 	1 -> wall point      -> [2,160]
! 		 	2 -> interior points -> {1},[161,9440]
! 		 	3 -> outer points    -> [9441,9600]

module data_structure_mod
!
!
	 use parameter_mod

	implicit none

	integer :: max_points,local_points
        integer::wall_points,interior_points,outer_points
		
	type :: points

!	!	scanned from input file	!	!
		real*8 :: x,y
		integer :: local_id
		integer :: global_id
		integer :: flag_1 ! stores location of point
		integer :: flag_2 ! stores shape point belongs to 
		integer :: nbhs
		integer, dimension(:),allocatable :: conn
!	!	!	!	!	!	!	!	!	!		

!		real*8 :: nx, ny
!		
!		real*8 :: rho, u1, u2, pr
!		real*8 :: flux_res(4)
!
!		real*8 :: q(4), qx(4), qy(4)
!
!		real*8 :: entropy, vorticity, vorticity_sqr

!		integer :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
!		integer :: xpos_conn(15), xneg_conn(15)
!		integer :: ypos_conn(15), yneg_conn(15)
!
	 end type points

         type(points), dimension(:), allocatable :: point


!
!
	
	save
!
	! integer :: wall_points_index(wall_points)
	! integer :: outer_points_index(outer_points)
	! integer :: interior_points_index(interior_points)	
	! integer :: shape_points_index(shapes, max_shape_points)
!
!

!       real*8	:: res_old, res_new, residue, max_res
!	integer :: max_res_point
!	real*8 	:: cfv
!	real*8	:: Cl, Cd, Cm
!	real*8	:: total_entropy, total_enstrophy


!
!
	end module data_structure_mod		 	 
