module data_structure_mod_diff
        
        use parameter_mod

        implicit none

        integer :: max_points,local_points,ghost_points
        integer :: wall_points,interior_points,outer_points,shape_points

!       ghost global indices
        integer , dimension(:), allocatable :: pghost

        type :: points

                integer, dimension(:), allocatable :: original_id
		real*8, dimension(:), allocatable :: x,y
                integer, dimension(:), allocatable :: left,right
                integer, dimension(:), allocatable :: flag_1 ! stores location of point
                integer, dimension(:), allocatable :: flag_2 ! stores shape point belongs to 
		real*8, dimension(:), allocatable :: nx,ny
                integer, dimension(:), allocatable :: nbhs
                integer, dimension(:,:), allocatable :: conn

                integer, dimension(:), allocatable :: qtdepth

                real*8, dimension(:), allocatable :: min_dist

                real*8, dimension(:,:), allocatable :: prim
                real*8, dimension(:,:), allocatable :: prim_old
		real*8, dimension(:,:), allocatable :: flux_res

                real*8, dimension(:,:), allocatable :: q
		real*8, dimension(:,:,:), allocatable :: dq
		real*8, dimension(:,:,:), allocatable :: qm

		real*8, dimension(:), allocatable :: entropy, vorticity, vorticity_sqr

                integer, dimension(:), allocatable :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer, dimension(:,:), allocatable :: xpos_conn, xneg_conn
                integer, dimension(:,:), allocatable :: ypos_conn, yneg_conn

                real*8, dimension(:), allocatable  :: delta
                
! Implicit data
                real*8, dimension(:,:), allocatable :: U_old
        end type points
 
        type :: points_diff

		real*8, dimension(:), allocatable :: nx,ny
		real*8, dimension(:), allocatable :: x,y
		real*8, dimension(:,:), allocatable :: x_res

                real*8, dimension(:,:), allocatable :: prim
                real*8, dimension(:,:), allocatable :: prim_old
		real*8, dimension(:,:), allocatable :: flux_res

                real*8, dimension(:,:), allocatable :: q
		real*8, dimension(:,:,:), allocatable :: dq
		real*8, dimension(:,:,:), allocatable :: qm

		real*8, dimension(:), allocatable :: entropy, vorticity, vorticity_sqr

                real*8, dimension(:), allocatable  :: delta
        end type points_diff
        type(points) :: point
        type(points_diff) :: pointb

        save

        integer,allocatable,dimension(:) :: wall_points_index
        integer,allocatable,dimension(:) :: outer_points_index
        integer,allocatable,dimension(:) :: interior_points_index
        integer,allocatable,dimension(:) :: shape_points_index

        !iterations
        integer :: it, itr, restart


        real*8  :: res_old, res_new, residue, max_res
        real* 8 :: gsum_res_sqr,sum_res_sqr
        integer :: max_res_point
	real*8, allocatable, dimension(:)  :: Cl, Cd, Cm, ClCd, cfv
	real*8, allocatable, dimension(:)  :: Clb, Cdb, Cmb, ClCdb
	real*8  :: total_entropy, total_enstrophy
	real*8  :: total_entropyb, total_enstrophyb
        integer :: plen

!The parameter CFL is the CFL number for stability ..
        real*8 :: CFL

        integer :: max_iters

! File format
        integer :: format
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
        integer :: obj_flag

!       No of shapes
        integer :: shapes

!       Adoint mode
        integer :: ad_mode

!       checkpoints
        integer :: chkpts

!       adjoint residue
        real*8 :: adj_res, adj_res_old

   contains

        subroutine allocate_soln()
                implicit none

                allocate(point%prim(4,max_points))
                allocate(point%prim_old(4,max_points))

                allocate(point%flux_res(4,max_points))

                allocate(point%U_old(4,max_points))

                allocate(point%q(4,max_points))

                allocate(point%dq(2,4,max_points))

                allocate(point%qm(2,4,max_points))

                allocate(point%entropy(max_points))
                allocate(point%vorticity(max_points))
                allocate(point%vorticity_sqr(max_points))

                
                allocate(point%xpos_nbhs(max_points))
                allocate(point%xneg_nbhs(max_points))
                allocate(point%ypos_nbhs(max_points))
                allocate(point%yneg_nbhs(max_points))


                allocate(point%xpos_conn(max_points,15))
                allocate(point%xneg_conn(max_points,15))


                allocate(point%ypos_conn(max_points,15))
                allocate(point%yneg_conn(max_points,15))


                allocate(point%delta(max_points))

                allocate(Cl(shapes))
                allocate(ClCd(shapes))
                allocate(Cd(shapes))
                allocate(Cm(shapes))
                allocate(cfv(shapes))

        end subroutine
        
        subroutine allocate_soln_b()
                implicit none

                allocate(pointb%prim(4,max_points))
                allocate(pointb%prim_old(4,max_points))

                allocate(pointb%flux_res(4,max_points))

                allocate(pointb%q(4,max_points))

                allocate(pointb%dq(2,4,max_points))

                allocate(pointb%qm(2,4,max_points))

                allocate(pointb%delta(max_points))

                allocate(Clb(shapes))
                allocate(ClCdb(shapes))
                allocate(Cdb(shapes))
                allocate(Cmb(shapes))

                allocate(pointb%vorticity_sqr(max_points))
        end subroutine

        subroutine deallocate_soln()
                implicit none

                deallocate(point%prim)
                deallocate(point%prim_old)

                deallocate(point%flux_res)
                deallocate(point%U_old)


                deallocate(point%q)
                deallocate(point%dq)
                deallocate(point%qm)


                deallocate(point%entropy)
                deallocate(point%vorticity)
                deallocate(point%vorticity_sqr)

                
                deallocate(point%xpos_nbhs)
                deallocate(point%xneg_nbhs)
                deallocate(point%ypos_nbhs)
                deallocate(point%yneg_nbhs)


                deallocate(point%xpos_conn)
                deallocate(point%xneg_conn)


                deallocate(point%ypos_conn)
                deallocate(point%yneg_conn)


                deallocate(point%delta)

                deallocate(Cl)
                deallocate(ClCd)
                deallocate(Cd)
                deallocate(Cm)
                deallocate(cfv)

        end subroutine

        subroutine deallocate_soln_b()
                implicit none

                deallocate(pointb%prim)
                deallocate(pointb%prim_old)

                deallocate(pointb%flux_res)

                deallocate(pointb%q)

                deallocate(pointb%dq)

                deallocate(pointb%qm)

                deallocate(pointb%delta)

                deallocate(Clb)
                deallocate(ClCdb)
                deallocate(Cdb)
                deallocate(Cmb)

                deallocate(pointb%vorticity_sqr)
        end subroutine
end module data_structure_mod_diff
