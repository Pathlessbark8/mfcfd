! flag_1 : 1 -> wall point      -> [2,160]
!          2 -> interior points -> {1},[161,9440]
! 	   3 -> outer points    -> [9441,9600]

module data_structure_mod
#include <petsc/finclude/petscvec.h>


        use petscvec
        use parameter_mod

        implicit none

        integer :: max_points,local_points,ghost_points
        integer :: wall_points,interior_points,outer_points
        integer :: shape_points(shapes)

!       ghost global indices
        integer , dimension(:), allocatable :: pghost


        type :: points

!	!	scanned from input file	!	!
		real*8, dimension(:), allocatable :: x,y
                integer, dimension(:), allocatable :: local_id
                integer, dimension(:), allocatable :: global_id
                integer, dimension(:), allocatable :: flag_1 ! stores location of point
                integer, dimension(:), allocatable :: flag_2 ! stores shape point belongs to 
                integer, dimension(:), allocatable :: nbhs
                integer, dimension(:,:), allocatable :: conn
!	!	!	!	!	!	!	

		real*8, dimension(:), allocatable :: nx, ny

		real*8, dimension(:,:), allocatable :: prim
		real*8, dimension(:,:), allocatable :: flux_res

		real*8, dimension(:,:), allocatable :: q
		real*8, dimension(:,:,:), allocatable :: dq

		real*8, dimension(:), allocatable :: entropy, vorticity, vorticity_sqr

                integer, dimension(:), allocatable :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer, dimension(:,:), allocatable :: xpos_conn, xneg_conn
                integer, dimension(:,:), allocatable :: ypos_conn, yneg_conn

                real*8, dimension(:), allocatable  :: delta

        end type points
 
        type(points) :: p

        save

        integer,allocatable,dimension(:) :: wall_points_index
        integer,allocatable,dimension(:) :: outer_points_index
        integer,allocatable,dimension(:) :: interior_points_index
!       TODO make below array dynamic for second index	
        integer :: shape_points_index(shapes, max_shape_points)

        !iterations
        integer*8 :: it


        real*8  :: res_old, res_new, residue, max_res
        integer :: max_res_point
	real*8  :: cfv
!	real*8	:: Cl, Cd, Cm
!	real*8	:: total_entropy, total_enstrophy
	real*8 :: total_enstrophy

!   PETSc variables

        PetscMPIInt              :: rank,proc
        Vec                  :: p_dq
        Vec                  :: p_prim

    contains

        subroutine allocate_soln()
                implicit none

                allocate(p%prim(4,max_points))


                allocate(p%flux_res(4,max_points))


                allocate(p%q(4,max_points))

                allocate(p%dq(2,4,max_points))


                allocate(p%entropy(max_points))
                allocate(p%vorticity(max_points))
                allocate(p%vorticity_sqr(max_points))

                
                allocate(p%xpos_nbhs(max_points))
                allocate(p%xneg_nbhs(max_points))
                allocate(p%ypos_nbhs(max_points))
                allocate(p%yneg_nbhs(max_points))


                allocate(p%xpos_conn(max_points,15))
                allocate(p%xneg_conn(max_points,15))


                allocate(p%ypos_conn(max_points,15))
                allocate(p%yneg_conn(max_points,15))


                allocate(p%delta(max_points))
        end subroutine


        subroutine init_petsc()
                implicit none
                PetscErrorCode       :: ierr
                if (proc==1) return
                if (rank==0) print*,'Setting up parallel vectors'
                pghost = pghost - 1

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,PETSC_DECIDE,ghost_points,pghost,p%dq(1,1,1),p_dq,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,PETSC_DECIDE,ghost_points,pghost,p%prim(1,1),p_prim,ierr)


        if(rank==0) print*,'Set up parallel vectors'
        end subroutine 

        subroutine dest_petsc()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
                    

                call VecDestroy(p_dq,ierr)
                call VecDestroy(p_prim,ierr)

        end subroutine 


        subroutine update_begin_dq_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine 

        subroutine update_begin_prim_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine 


        subroutine update_end_dq_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine 

        subroutine update_end_prim_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine


end module data_structure_mod
