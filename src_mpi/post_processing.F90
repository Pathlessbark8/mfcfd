module post_processing_mod
	

        use data_structure_mod
	

contains


	subroutine print_primal_output()
	
	
		implicit none
			
		integer :: k
		
		OPEN(UNIT=101,FILE="primal-solution.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		
		
		write(101,*) 'TITLE="Transonic flow over NACA0012 airfoil"'
		write(101,*) 'Variables="x","y","rho","u1","u2","pr","sensor"'
		write(101,*) 'Zone I = 160, J = 60, F = POINT'
			do k=1, max_points
				write(101,*) point%x(k), point%y(k), point%prim(4,k), point%prim(2,k), point%prim(3,k), &
				point%prim(1,k),  point%sensor(k)
			enddo
					
		CLOSE(UNIT=101)
			
	end subroutine		




end module post_processing_mod
