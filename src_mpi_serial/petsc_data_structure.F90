module petsc_data_structure_mod
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petsclog.h>
    
    
    use petscvec
    use parameter_mod
    use data_structure_mod
    
    implicit none
    
    PetscMPIInt          :: rank,proc
    Vec                  :: p_dq
    Vec                  :: p_ddq
    Vec                  :: p_qm
    Vec                  :: p_prim
    Vec                  :: p_du1_dx
    Vec                  :: p_du2_dy
    PetscLogEvent        :: dq_comm, ddq_comm, prim_comm, qm_comm, du1_dx_comm, du2_dy_comm
    
    contains
    
    subroutine init_petsc()
        implicit none
        PetscErrorCode       :: ierr
        if(rank == 0) then
            write(*,*) '%%%%%%%-Setting up parallel vectors-%%%%%%%'
            write(*,*)
        end if
        pghost = pghost - 1
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%dq(1,1,1),p_dq,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%qm(1,1,1),p_qm,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,3*4,3*4*local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%ddq(1,1,1),p_ddq,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%prim(1,1),p_prim,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%du1_dx(1),p_du1_dx,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
        &PETSC_DECIDE,ghost_points,pghost,point%du2_dy(1),p_du2_dy,ierr)
        
        call VecGetSize(p_prim,plen,ierr)
        plen = plen/4
        
        call PetscLogEventRegister('dq_comm',  0,dq_comm,ierr);
        call PetscLogEventRegister('ddq_comm',  0,ddq_comm,ierr);
        call PetscLogEventRegister('qm_comm',  0,qm_comm,ierr);
        call PetscLogEventRegister('prim_comm',  0,prim_comm,ierr);
        call PetscLogEventRegister('du1_dx_comm',  0,du1_dx_comm,ierr);
        call PetscLogEventRegister('du2_dy_comm',  0,du2_dy_comm,ierr);
        
    end subroutine 
    
    subroutine dest_petsc()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        
        call VecDestroy(p_dq,ierr)
        call VecDestroy(p_ddq,ierr)
        call VecDestroy(p_qm,ierr)
        call VecDestroy(p_prim,ierr)
        call VecDestroy(p_du1_dx,ierr)
        call VecDestroy(p_du2_dy,ierr)
        
    end subroutine 
    
    
    subroutine update_begin_dq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call VecGhostUpdateBegin(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        
    end subroutine 
    
    subroutine update_begin_ddq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call VecGhostUpdateBegin(p_ddq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        
    end subroutine 
    
    subroutine update_begin_prim_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call VecGhostUpdateBegin(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
        
    end subroutine 
    
    subroutine update_begin_qm_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call VecGhostUpdateBegin(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
        
    end subroutine

    subroutine update_begin_du1_dx_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        call VecGhostUpdateBegin(p_du1_dx,INSERT_VALUES,SCATTER_FORWARD,ierr)
    end subroutine 

    subroutine update_begin_du2_dy_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        call VecGhostUpdateBegin(p_du2_dy,INSERT_VALUES,SCATTER_FORWARD,ierr)
    end subroutine 
    
    subroutine update_end_dq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(dq_comm, ierr)
        call VecGhostUpdateEnd(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(dq_comm, ierr)
        
    end subroutine 
    
    subroutine update_end_ddq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(ddq_comm, ierr)
        call VecGhostUpdateEnd(p_ddq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(ddq_comm, ierr)
        
    end subroutine 
    
    subroutine update_end_qm_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(qm_comm, ierr)
        call VecGhostUpdateEnd(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(qm_comm, ierr)
        
    end subroutine 
    
    subroutine update_end_prim_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(prim_comm, ierr)
        call VecGhostUpdateEnd(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(prim_comm, ierr)
        
    end subroutine
    
    subroutine update_end_du1_dx_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(du1_dx_comm, ierr)
        call VecGhostUpdateEnd(p_du1_dx,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(du1_dx_comm, ierr)
        
    end subroutine

    subroutine update_end_du2_dy_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
        
        call PetscLogEventBegin(du2_dy_comm, ierr)
        call VecGhostUpdateEnd(p_du2_dy,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(du2_dy_comm, ierr)
        
    end subroutine
    
end module petsc_data_structure_mod
