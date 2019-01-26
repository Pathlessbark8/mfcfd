module fpi_solver_mod
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

                call eval_q_variables()
                
                call eval_q_derivatives()

                call func_delta()   

                call cal_flux_residual()

                call state_update()

                
                call objective_function()

                res_new = dsqrt(sum_res_sqr)/plen

                if(t .le. 2) then
                        res_old = res_new
                        residue = 0.d0
                else 
                        residue = dlog10(res_new/res_old)
                endif

        end subroutine

end module fpi_solver_mod
