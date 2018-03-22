module q_lskum_mod

!	First written on 14.10.2016
!	updated on Dec 26, 2016
!	updated on Dec 29, 2016

	use data_structure_mod
	use point_normals_mod    
	use generate_connectivity_mod
	use fpi_solver_mod	

contains

	subroutine q_lskum()

		implicit none
	
		integer :: t, i
		external cpu_time
		real*8 :: start_time, end_time
		OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!		call cpu_time(start_time)

		call compute_normals()
		call generate_connectivity()

			
		do t = 1, 5
			call fpi_solver(t)
!			print*, t, res_new, residue, max_res, max_res_point, Cl, Cd, cfv
!			print*, t, res_new, residue, Cl, Cd, cfv
!			write(301, *) t, residue				
		enddo

!		call cpu_time(end_time)
!		print*, "net time = ", end_time - start_time

		CLOSE(UNIT=301)	
							
	end subroutine

end module q_lskum_mod
