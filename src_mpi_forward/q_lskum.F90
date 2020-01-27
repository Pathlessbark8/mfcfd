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

                integer :: i

                OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

                call compute_normals()
                call generate_connectivity()

                write(*,*)'%%%%-Normals and connectivity generated-%%%'
                write(*,*)

                do i=1,max_points
                        point%phi1(:,i) = 1.0d0
                        point%phi2(:,i) = 1.0d0
                enddo

                ! point%phi1(80,1) = point%phi1(80,1) + 1e-3
        
                write(*,*)'%%%%%%%%%%%%%-Iterations begin-%%%%%%%%%%%%'
                write(*,*)

                if(restart == 0)itr = 0
                
                do it = itr+1, itr+max_iters
                        
                        call fpi_solver(it)

                        write(*,'(a12,i8,a15,e30.20)')'iterations:',it,'residue:',residue
                        write(301, *) it, residue
                        if(residue.ne.residue)exit
                
                enddo
                
                CLOSE(UNIT=301)

        end subroutine

end module q_lskum_mod
