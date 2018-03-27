! flag_1 : 	1 -> wall point      -> [2,160]
! 		 	2 -> interior points -> {1},[161,9440]
! 		 	3 -> outer points    -> [9441,9600]

module data_structure_mod
#include <petsc/finclude/petscvec.h>


        use petscvec
        use parameter_mod

        implicit none

        integer :: max_points,local_points,ghost_points
        integer :: wall_points,interior_points,outer_points
        integer :: shape_points(shapes)

!   ghost global indices
        integer , dimension(:), allocatable :: pghost

! !   data structure to hold points by location	
! 	integer , dimension(:), allocatable :: wall_points_index
! 	integer , dimension(:), allocatable :: interior_points_index
! 	integer , dimension(:), allocatable :: outer_points_index

! !	data structure to hold points by shape
! 	integer, dimension(:,:),allocatable :: shape_points_index


        type :: points

!	!	scanned from input file	!	!
		real*8 :: x,y
                integer :: local_id
                integer :: global_id
                integer :: flag_1 ! stores location of point
                integer :: flag_2 ! stores shape point belongs to 
                integer :: nbhs
                integer :: conn(15)
!	!	!	!	!	!	!	!	!	!		

		real*8 :: nx, ny

		real*8 :: rho, u1, u2, pr
		real*8 :: flux_res(4)

		real*8 :: q(4), qx(4), qy(4)

		real*8 :: entropy, vorticity, vorticity_sqr

                integer :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer :: xpos_conn(15), xneg_conn(15)
                integer :: ypos_conn(15), yneg_conn(15)

        end type points
 
        type(points), dimension(:), allocatable :: point

        save

        integer,allocatable,dimension(:) :: wall_points_index
        integer,allocatable,dimension(:) :: outer_points_index
        integer,allocatable,dimension(:) :: interior_points_index
!TODO make below array dynamic for second index	
        integer :: shape_points_index(shapes, max_shape_points)


!   real*8	:: res_old, res_new, residue, max_res
!	integer :: max_res_point
!	real*8 	:: cfv
!	real*8	:: Cl, Cd, Cm
!	real*8	:: total_entropy, total_enstrophy

!   PETSc variables

        integer              :: rank,proc
        Vec                  :: p_q1,p_q2,p_q3,p_q4
        Vec                  :: px_q1,px_q2,px_q3,px_q4
        Vec                  :: py_q1,py_q2,py_q3,py_q4
        Vec                  :: p_u1, p_u2, p_pr, p_rho

    contains

        subroutine init_petsc()
                implicit none
                PetscErrorCode       :: ierr
                if (proc==1) return
                if (rank==0) print*,'Setting up parallel vectors'
                ghost_points = ghost_points - 1

                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%q(1),p_q1,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%q(2),p_q2,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%q(3),p_q3,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%q(4),p_q4,ierr)

                
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qx(1),px_q1,ierr) 
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qx(2),px_q2,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qx(3),px_q3,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qx(4),px_q4,ierr)


                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qy(1),py_q1,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qy(2),py_q2,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qy(3),py_q3,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%qy(4),py_q4,ierr)
                

                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%u1, p_u1,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%u2, p_u2,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%pr, p_pr,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,point(:)%rho,p_rho,ierr)
        

        if(rank==0) print*,'Set up parallel vectors'
        end subroutine init_petsc

        subroutine dest_petsc()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
                    

                call VecDestroy(p_q1,ierr)
                call VecDestroy(p_q2,ierr)
                call VecDestroy(p_q3,ierr)
                call VecDestroy(p_q4,ierr)

        end subroutine dest_petsc

        subroutine update_begin_q_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_q_ghost


        subroutine update_begin_qx_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(px_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(px_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(px_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(px_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_qx_ghost


        subroutine update_begin_qy_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(py_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(py_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(py_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(py_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_qy_ghost

        subroutine update_begin_u1_u2_pr_rho_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_u1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_u2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_pr,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(p_rho,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_u1_u2_pr_rho_ghost

        subroutine update_end_q_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_q_ghost


        subroutine update_end_qx_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(px_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(px_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(px_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(px_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_qx_ghost


        subroutine update_end_qy_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(py_q1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(py_q2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(py_q3,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(py_q4,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_qy_ghost


        subroutine update_end_u1_u2_pr_rho_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_u1,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_u2,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_pr,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(p_rho,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_u1_u2_pr_rho_ghost







end module data_structure_mod
