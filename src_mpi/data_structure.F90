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

		real*8, dimension(:), allocatable :: rho, u1, u2, pr
		real*8, dimension(:,:), allocatable :: flux_res

		real*8, dimension(:,:), allocatable :: q, qx, qy

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


        real*8  :: res_old, res_new, residue, max_res
        integer :: max_res_point
!	real*8 	:: cfv
!	real*8	:: Cl, Cd, Cm
!	real*8	:: total_entropy, total_enstrophy

!   PETSc variables

        PetscMPIInt              :: rank,proc
        Vec                  :: p_q
        Vec                  :: p_qx
        Vec                  :: p_qy
        Vec                  :: p_u1, p_u2, p_pr, p_rho

    contains

        subroutine allocate_soln()
                implicit none

                allocate(p%nx(max_points))
                allocate(p%ny(max_points))
                

                allocate(p%rho(max_points))
                allocate(p%u1(max_points))
                allocate(p%u2(max_points))
                allocate(p%pr(max_points))


                allocate(p%flux_res(4,max_points))


                allocate(p%q(4,max_points))
                allocate(p%qx(4,max_points))
                allocate(p%qy(4,max_points))


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

      
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,PETSC_DECIDE,ghost_points,pghost,p%q(1,1),p_q,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,PETSC_DECIDE,ghost_points,pghost,p%qx(1,1),p_qx,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,PETSC_DECIDE,ghost_points,pghost,p%qy(1,1),p_qy,ierr)

                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,p%u1(1), p_u1,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,p%u2(1), p_u2,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,p%pr(1), p_pr,ierr)
                call VecCreateGhostWithArray(PETSC_COMM_WORLD,local_points,PETSC_DECIDE,ghost_points,pghost,p%rho(1),p_rho,ierr)
        

        if(rank==0) print*,'Set up parallel vectors'
        end subroutine init_petsc

        subroutine dest_petsc()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
                    

                call VecDestroy(p_q,ierr)
                call VecDestroy(p_qx,ierr)
                call VecDestroy(p_qy,ierr)
                call VecDestroy(p_u1,ierr)
                call VecDestroy(p_u2,ierr)
                call VecDestroy(p_pr,ierr)
                call VecDestroy(p_rho,ierr)

        end subroutine dest_petsc

        subroutine update_begin_q_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_q,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_q_ghost


        subroutine update_begin_qx_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_qx,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_begin_qx_ghost


        subroutine update_begin_qy_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_qy,INSERT_VALUES,SCATTER_FORWARD,ierr)

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

                call VecGhostUpdateEnd(p_q,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_q_ghost


        subroutine update_end_qx_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_qx,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine update_end_qx_ghost


        subroutine update_end_qy_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(p_qy,INSERT_VALUES,SCATTER_FORWARD,ierr)

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
