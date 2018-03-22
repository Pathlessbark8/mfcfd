module fpi_solver_mod
#include <petsc/finclude/petscsys.h>
!
!	First written on 14.10.2016
!	updated on Dec 26, 2016
!	updated on Dec 29, 2016
!	updated on Oct 30, 2017	
!	
!
	use data_structure_mod
!	use flux_residual_mod
!	use state_update_mod
	use q_variables_mod
!	use objective_function_mod
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
                        !Update the ghost values from the owned process
                        call update_begin_ghost()

			call eval_q_variables()

                        !End the update of ghost values
                        call update_end_ghost()


			call eval_q_derivatives()	

!				
!			call cal_flux_residual()
!			call state_update()
!			call objective_function()
!			
!			if(t .le. 2) then
!					res_old = res_new
!					residue = 0.d0
!			else 
!					residue = dlog10(res_new/res_old)
!			endif					
!							
	end subroutine
!
!	
!
end module fpi_solver_mod
