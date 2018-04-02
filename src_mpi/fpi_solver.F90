module fpi_solver_mod
#include <petsc/finclude/petscsys.h>
!
!	First written on 14.10.2016
!	updated on Dec 26, 2016
!	updated on Dec 29, 2016
!	updated on Oct 30, 2017	
!	
!
        use data_structure_mod
        use flux_residual_mod
        use state_update_mod
        use q_variables_mod
        use objective_function_mod


contains


        subroutine fpi_solver(t)

                implicit none
                
                integer*8 :: t
                real*8 :: lresidue
                PetscErrorCode :: ierr

                call eval_q_variables()
                
                call eval_q_derivatives()

                !Update the ghost values from the owned process

                call update_begin_dq_ghost()

                call update_end_dq_ghost()


                call cal_flux_residual()

                call func_delta()            

                call state_update()

                call update_begin_prim_ghost()
                call update_end_prim_ghost()

                call objective_function()


                if(t .le. 2) then
                        res_old = res_new
                        lresidue = 0.d0
                else 
                        lresidue = dlog10(res_new/res_old)
                endif
                call MPI_Reduce(lresidue, residue, 1, MPI_DOUBLE, MPI_SUM, &
                   0, PETSC_COMM_WORLD, ierr)

        
        end subroutine

end module fpi_solver_mod
