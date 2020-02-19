module objective_function_mod
    
    use data_structure_mod
    use compute_force_coeffs_mod
    use compute_entropy_mod
    use compute_enstrophy_mod
    use stagnation_values_mod

    contains

        subroutine objective_function()

            implicit none

            call compute_cl_cd_cm()
            call compute_entropy()
            call compute_enstrophy()
            call objective_function_J()
            vector_cost_func = total_loss_stagpressure + Cl + Cd + Cm + ClCd + total_entropy + total_enstrophy
            if(rank==0) then
                write(*,*) "SG", total_loss_stagpressure 
                write(*,*) "Cl", cl
                write(*,*) "Cd", cd
                write(*,*) "Cm", cm
                write(*,*) "Clcd", clcd
                write(*,*) "Entropy", total_entropy
                write(*,*) "Enstrophy", total_enstrophy
                write(*,*) "Vector function is ", vector_cost_func
              end if
        end subroutine objective_function

end module objective_function_mod
