module petsc_data_structure_mod
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petsclog.h>


        use petscvec
        use parameter_mod
        use data_structure_mod_diff

        implicit none

        PetscMPIInt          :: rank,proc
        Vec                  :: p_dq
        Vec                  :: p_qm
        Vec                  :: p_prim
        Vec                  :: pd_dq
        Vec                  :: pd_qm
        Vec                  :: pd_prim

        ! Adjoint PETSc variables
        Vec                  :: pb_dq, pb_q, pb_qm
        Vec                  :: pb_prim
        Vec                  :: pdb_dq, pdb_q, pdb_qm
        Vec                  :: pdb_prim
        Vec                  :: pb_x, pb_y
        PetscLogEvent        :: dq_comm, prim_comm, qm_comm

    contains

        subroutine init_petsc()
                implicit none
                PetscErrorCode       :: ierr
                if(rank == 0) then
                        write(*,*) '%%%%%%%-Setting up parallel vectors-%%%%%%%'
                        write(*,*)
                end if
                pghost = pghost - 1

                ! Primal variables
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,point%dq(1,1,1),p_dq,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,point%qm(1,1,1),p_qm,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,point%prim(1,1),p_prim,ierr)
                
                ! Adjoint variables
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%dq(1,1,1),pb_dq,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%prim(1,1),pb_prim,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%q(1,1),pb_q,ierr)

                call VecCreateGhostWithArray(PETSC_COMM_WORLD, local_points, &
                                 &PETSC_DECIDE, ghost_points, pghost, pointb%x(1), &
                                 &pb_x, ierr); CHKERRQ(ierr)

                call VecCreateGhostWithArray(PETSC_COMM_WORLD, local_points, &
                                 &PETSC_DECIDE, ghost_points, pghost, pointb%y(1), &
                                 &pb_y, ierr); CHKERRQ(ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%qm(1,1,1),pb_qm,ierr)

                ! forward part in the adjoint code
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointd%dq(1,1,1),pd_dq,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointd%qm(1,1,1),pd_qm,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointd%prim(1,1),pd_prim,ierr)

                ! forward adjoint part in the aoa adjoint code
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointdb%dq(1,1,1),pdb_dq,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointdb%prim(1,1),pdb_prim,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointdb%q(1,1),pdb_q,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointdb%qm(1,1,1),pdb_qm,ierr)
                
                call VecGetSize(p_prim,plen,ierr)
                plen = plen/4

                call PetscLogEventRegister('dq_comm',  0,dq_comm,ierr);
                call PetscLogEventRegister('qm_comm',  0,qm_comm,ierr);
                call PetscLogEventRegister('prim_comm',  0,prim_comm,ierr);

        end subroutine 

        subroutine dest_petsc()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecDestroy(p_dq,ierr)
                call VecDestroy(p_qm,ierr)
                call VecDestroy(p_prim,ierr)
                call VecDestroy(pb_dq,ierr)
                call VecDestroy(pb_prim,ierr)
                call VecDestroy(pb_q,ierr)
                call VecDestroy(pb_x,ierr)
                call VecDestroy(pb_y,ierr)

                call VecDestroy(pd_dq,ierr)
                call VecDestroy(pd_qm,ierr)
                call VecDestroy(pd_prim,ierr)
                
                call VecDestroy(pdb_dq,ierr)
                call VecDestroy(pdb_qm,ierr)
                call VecDestroy(pdb_prim,ierr)
                call VecDestroy(pdb_q,ierr)
        end subroutine 

        subroutine update_begin_dq_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(pd_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine

        subroutine update_begin_prim_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(pd_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine

        subroutine update_begin_qm_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateBegin(pd_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)

        end subroutine

        subroutine update_end_dq_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(dq_comm, ierr)
                call VecGhostUpdateEnd(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(pd_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventEnd(dq_comm, ierr)

        end subroutine

        subroutine update_end_qm_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(qm_comm, ierr)
                call VecGhostUpdateEnd(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(pd_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventEnd(qm_comm, ierr)

        end subroutine

        subroutine update_end_prim_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(prim_comm, ierr)
                call VecGhostUpdateEnd(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call VecGhostUpdateEnd(pd_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventEnd(prim_comm, ierr)

        end subroutine

        subroutine update_begin_dqb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(pb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateBegin(pdb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 
        
        subroutine update_begin_primb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateBegin(pdb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 

        subroutine update_begin_qb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_q,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateBegin(pdb_q,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 
        
        subroutine update_begin_xb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_x,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 
       
        subroutine update_begin_yb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_y,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 
        
        subroutine update_begin_qmb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(pb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateBegin(pdb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 
        
        subroutine update_end_dqb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateEnd(pdb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 

        subroutine update_end_primb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateEnd(pdb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine
        
        subroutine update_end_qb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_q,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateEnd(pdb_q,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine
        
        subroutine update_end_xb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_x,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine
        
        subroutine update_end_yb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_y,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine
        
        subroutine update_end_qmb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateEnd(pb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)
                call VecGhostUpdateEnd(pdb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)

        end subroutine 

end module petsc_data_structure_mod
