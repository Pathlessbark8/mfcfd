module fpi_solver_mod
        use data_structure_mod
        use flux_residual_mod
        use state_update_mod
        use q_variables_mod
        use objective_function_mod


contains


        subroutine fpi_solver(t)

                implicit none
                
                integer :: t, i, rk

                do i =1, local_points
                        point%prim_old(:, i) = point%prim(:, i)
                end do

                call eval_q_variables()
                
                call eval_q_derivatives()

                call func_delta()   
                
                ! Perform 4-stage, 3-order SSPRK update
                do rk = 1, 4

                        call cal_flux_residual()

                        call state_update(rk)
                        
                end do

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
