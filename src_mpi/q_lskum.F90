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

                OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!TODO               Temporarily, until an effiecient way of computing normals in parallel is developed
!                call compute_normals()
                call generate_connectivity()

        
                do it = 1, max_iters
                        
                        call fpi_solver(it)
!                        if (rank==0) write(*,'(a12,i8,a15,e30.20)')'iterations:',it,'residue:',residue
!			print*, t, res_new, residue, max_res, max_res_point, Cl, Cd, cfv
                        write(301, *) it, residue
                enddo

!		print*, "net time = ", end_time - start_time

!		CLOSE(UNIT=301)	

        end subroutine

end module q_lskum_mod
