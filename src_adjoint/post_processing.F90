module post_processing_mod

        use data_structure_mod_diff
        use petsc_data_structure_mod

contains

        subroutine print_output()

                implicit none

                integer :: i
                character(len=64) :: sfile
                character(len=10) :: itos

                if(proc==1) then
                        sfile = 'solution/sol.dat'
                else
                        sfile = 'solution/'//'sol-'//trim(itos(4,rank))//'.dat'
                end if

                OPEN(UNIT=501,FILE=trim(sfile))
                
                write(501,*)local_points

                do i = 1, local_points
                        write(501,'(i6, 4e30.20)')point%original_id(i), point%x(i),&
                                point%y(i),pointb%x(i), pointb%y(i)
                end do

                close(501)



        end subroutine

end module post_processing_mod
