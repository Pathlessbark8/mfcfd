module initial_conditions_mod
!
!
	use data_structure_mod
!
!
contains
!
!
	subroutine initial_conditions()	
	!
	!
		implicit none
	!
		
		integer :: k,i
		!
		!
		if(initial_conditions_flag .eq. 0) then
		!
				do k=1, local_points
		!
					point(k)%rho =	rho_inf
					point(k)%u1	= Mach*dcos(theta)
					point(k)%u2	= Mach*dsin(theta)
					point(k)%pr	= pr_inf
				enddo		
!		!	
		else if(initial_conditions_flag .eq. 1) then
		!
				OPEN(UNIT=100,FILE="./structured-grid-160-60/stored-solution",FORM="FORMATTED",STATUS="OLD",ACTION="READ")

				do k=1, local_points
						read(100,*) i, point(k)%rho, point(k)%u1, point(k)%u2, point(k)%pr
				enddo	
	!
				CLOSE(UNIT=100)
		endif		
	!
	end subroutine
!
!
end module initial_conditions_mod
!		

