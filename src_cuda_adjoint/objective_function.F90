module objective_function_mod
        
        use DATA_STRUCTURE_MOD_DIFF
        use compute_force_coeffs_mod
        use compute_entropy_mod

        contains


                subroutine objective_function()


                        implicit none
                        
                        call compute_cl_cd_cm()
                        call compute_entropy()

                end subroutine objective_function

end module objective_function_mod
