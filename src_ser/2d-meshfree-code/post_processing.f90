module post_processing_mod
!
!	First written on 14.10.2016
!	updated on April 15, 2017
!
	
!
	use data_structure_mod
!	
!
contains
!
!
	subroutine print_primal_output()
	!
	!
		implicit none
!			
		integer :: k
		
		OPEN(UNIT=101,FILE="primal-solution.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		
		
		write(101,*) 'TITLE="Transonic flow over NACA0012 airfoil"'
		write(101,*) 'Variables="x","y","rho","u1","u2","pr","entropy","vorticity","vorticity_sqr"'
		write(101,*) 'Zone I = 160, J = 60, F = POINT'
			do k=1, max_points
				write(101,*) point(k)%x, point(k)%y, point(k)%rho, point(k)%u1, point(k)%u2, &
				point(k)%pr, point(k)%entropy, point(k)%vorticity, point(k)%vorticity_sqr
			enddo
	!				
		CLOSE(UNIT=101)
	!		
	end subroutine		
!
!
!
!
end module post_processing_mod
