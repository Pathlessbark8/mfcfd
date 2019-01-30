module objective_function_mod
        
        use data_structure_mod
        use compute_force_coeffs_mod

        contains


                subroutine objective_function()


                        implicit none
                        
                        call compute_cl_cd_cm()

                end subroutine objective_function
!
!									
!
end module objective_function_mod
!		
!
!
