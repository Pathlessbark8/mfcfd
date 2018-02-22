module objective_function_mod
!
	use data_structure_mod
	use compute_force_coeffs_mod
	use compute_entropy_mod
	use compute_enstrophy_mod
!
	contains
!
!
		subroutine objective_function()
!
!
			implicit none
!
!
			call compute_cl_cd_cm()
			call compute_entropy()
			call compute_enstrophy()
!
!			CFV = Cl_flag*Cl + Cd_flag*Cd + Cm_flag*Cm + Cl_Cd_flag*Cl/Cd
!			CFV = Cl/Cd
!			CFV = total_entropy
			CFV = total_enstrophy
!
!
		end subroutine objective_function
!
!									
!
end module objective_function_mod
!		
!
!
