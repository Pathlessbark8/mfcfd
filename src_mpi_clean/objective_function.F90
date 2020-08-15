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
            ! call compute_entropy()
            ! call compute_enstrophy()
            ! call objective_function_J()
            vector_cost_func = Cl
        end subroutine objective_function

end module objective_function_mod
