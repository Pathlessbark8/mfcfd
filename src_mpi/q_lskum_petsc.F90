module q_lskum_petsc_mod
#include <petsc/finclude/petscts.h>

        use petscts
        use data_structure_mod
        use petsc_data_structure_mod
        use initial_conditions_mod
        use point_normals_mod    
        use generate_connectivity_mod
        use state_update_mod
        use q_variables_mod
        use flux_residual_mod

contains

        subroutine q_lskum_petsc()
!
!                implicit none
!
!                integer :: i
!
!                TS :: ts
!
!                PetscErrorCode :: ierr
!
!                if(rank==0)OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!
!!	        Assign the initial conditions for the primitive variables ..	
!
!                call initial_conditions()
!                if(rank == 0) then
!                        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
!                        write(*,*)
!                end if
!
!                call compute_normals()
!                call generate_connectivity()
!
!                if(rank == 0) then
!                        write(*,*)'%%%%-Normals and connectivity generated-%%%'
!                        write(*,*)
!                end if
!
!                call func_delta()
!
!                call TSCreate(PETSC_COMM_WORLD, ts, ierr); CHKERRQ(ierr)
!                
!                call TSSetProblemType(ts, TS_NONLINEAR, ierr); CHKERRQ(ierr)
!
!
!                
!                CLOSE(UNIT=301)
!
        end subroutine
!
!        subroutine fpi_solver()
!                
!                implicit none
!
!                call update_begin_prim_ghost()
!                call update_end_prim_ghost()
!
!                call prim_to_conserved()
!                        
!                call eval_q_variables()
!                
!                call eval_q_derivatives()
!                
!                !Update the ghost values from the owned process
!                call update_begin_dq_ghost()
!                call update_begin_qm_ghost()
!                call update_end_dq_ghost()
!                call update_end_qm_ghost()
!                
!                call cal_flux_residual()
!
!        end subroutine
!
!        subroutine prim_to_conserved()
!                
!                implicit none
!
!                integer :: i, k
!                real*8 :: nx, ny
!
!                do i = 1, wall_points
!                        
!                        k = wall_points_index(i)
!
!                        nx = point%nx(k)
!                        ny = point%ny(k)
!
!                        call primitive_to_conserved(point%prim(:,k), nx, ny, point%U(:,k))
!                end do
!
!                do i = 1, outer_points
!                        
!                        k = outer_points_index(i)
!
!                        nx = point%nx(k)
!                        ny = point%ny(k)
!                        
!                        call conserved_vector_Ubar(point%prim(:,k), point%U(:,k), nx, ny)
!                end do
!
!                do i = 1, interior_points
!                        
!                        k = interior_points_index(i)
!
!                        nx = point%nx(k)
!                        ny = point%ny(k)
!
!                        call primitive_to_conserved(point%prim(:,k), nx, ny, point%U(:,k))
!                end do
!
!
!        end subroutine

end module q_lskum_petsc_mod
