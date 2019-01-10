module fpi_solver_mod
#include <petsc/finclude/petscsys.h>
        use data_structure_mod
        use flux_residual_mod
        use state_update_mod
        use q_variables_mod
        use objective_function_mod
        use post_processing_mod 


contains


        subroutine fpi_solver(t)

                implicit none
                
                integer :: t, i
                PetscErrorCode :: ierr


                call eval_q_variables()
                
                call eval_q_derivatives()

                !Update the ghost values from the owned process

                call update_begin_dq_ghost()
                call func_delta()   
                call update_end_dq_ghost()

                call cal_flux_residual()

                ! Save previous solution
                do i=1,max_points
                        call primitive_to_conserved(i, point%nx(i), point%ny(i), point%U_old(:,i))
                end do
                
                call state_update()

                call update_begin_prim_ghost()
                call update_end_prim_ghost()

                call objective_function()

                call MPI_Reduce(sum_res_sqr,gsum_res_sqr, 1, MPI_DOUBLE, MPI_SUM, &
                   0, PETSC_COMM_WORLD, ierr)

                call MPI_Bcast(gsum_res_sqr, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD, &
                  ierr)

                res_new = dsqrt(gsum_res_sqr)/plen

                if(t .le. 2) then
                        res_old = res_new
                        residue = 0.d0
                else 
                        residue = dlog10(res_new/res_old)
                endif

                ! Primal output
                if(mod(it,nsave)==0)call print_primal_output()

        
        end subroutine

end module fpi_solver_mod
