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
                
                integer :: t, i, rk
                PetscErrorCode :: ierr

                do i =1, local_points
                        point%prim_old(:, i) = point%prim(:, i)
                end do

                call func_delta()   

                ! Perform 4-stage, 3-order SSPRK update
                do rk = 1, rks
                
                        call eval_q_variables()
                
                        call eval_q_derivatives()
                
                        !Update the ghost values from the owned process
                        call update_begin_dq_ghost()
                        call update_begin_qm_ghost()
                        call update_end_dq_ghost()
                        call update_end_qm_ghost()

                        call eval_q_double_derivatives()

                        call update_begin_ddq_ghost()
                        call update_end_ddq_ghost()

                        if(inner_iterations /= 0) then
                                do i = 1, inner_iterations
                                        call eval_dq_inner_loop()
                                        call eval_update_innerloop_2()
                                        call update_begin_dq_ghost()
                                        call update_end_dq_ghost()
                                        call eval_ddq_inner_loop()
                                        call eval_update_innerloop_3()
                                        call update_begin_ddq_ghost()
                                        call update_end_ddq_ghost()
                                enddo
                        end if
                
                        call cal_flux_residual()

                        call state_update(rk)

                        ! start updating primitive values
                        call update_begin_prim_ghost()
                        call update_end_prim_ghost()
                end do
                

                call objective_function()

                call MPI_Reduce(sum_res_sqr,gsum_res_sqr, 1, MPI_DOUBLE, MPI_SUM, &
                   0, PETSC_COMM_WORLD, ierr)

                call MPI_Bcast(gsum_res_sqr, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD, &
                  ierr)

                res_new = dsqrt(gsum_res_sqr)/plen

                if(t .le. 2 .and. restart == 0) then
                        res_old = res_new
                        residue = 0.d0
                else 
                        residue = dlog10(res_new/res_old)
                endif

                ! Print primal output
                if(mod(it,nsave)==0) then
                        if(rank == 0) then
                                write(*,*)
                                write(*,*)'%%%%%%%%%%%%%-Saving solution-%%%%%%%%%%%%%'
                                write(*,*)
                        end if
                        call print_primal_output()
                end if
                

        end subroutine

end module fpi_solver_mod
