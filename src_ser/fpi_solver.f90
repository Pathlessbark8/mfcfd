module fpi_solver_mod
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
!
!
contains
!
!
	subroutine fpi_solver(t)
!
		implicit none
		!		
			integer :: t,i 
!
			call eval_q_variables()
			call eval_q_derivatives()	

!				
			call cal_flux_residual()
			call state_update()
			call objective_function()
!			
			if(t .le. 2) then
					res_old = res_new
					residue = 0.d0
			else 
					residue = dlog10(res_new/res_old)
			endif					
!							
	end subroutine
!
!	
!
end module fpi_solver_mod
