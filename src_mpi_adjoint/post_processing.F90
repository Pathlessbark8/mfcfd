module post_processing_mod


        use data_structure_mod_diff
        use petsc_data_structure_mod


contains


        subroutine print_primal_output()


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
                
                write(501,*)local_points, it-1 , res_old

                do i = 1, local_points
                        write(501,'(3i8,6e30.20)')point%original_id(i), &
                                & point%flag_1(i),point%flag_2(i),point%x(i), &
                                & point%y(i),point%prim(1,i),point%prim(2,i),point%prim(3,i), &
                                & point%prim(4,i)
                end do

                close(501)

        end subroutine

        subroutine print_phi_output()


                implicit none

                integer :: i
                character(len=64) :: pfile
                character(len=10) :: itos

                if(proc==1) then
                        pfile = 'phistore/phid.dat'
                else
                        pfile = 'phistore/'//'phid-'//trim(itos(4,rank))//'.dat'
                end if

                OPEN(UNIT=501,FILE=trim(pfile))
                
                write(501,*)local_points, it-1 , res_old

                do i = 1, max_points
                        write(501,'(1i8,8e30.20)')point%original_id(i), &
                                & pointb%phi1(1,i), pointb%phi1(2,i), pointb%phi1(3,i), &
                                & pointb%phi1(4,i), pointb%phi2(1,i), pointb%phi2(2,i), pointb%phi2(3,i), &
                                & pointb%phi2(4,i)
                end do

                close(501)

        end subroutine

end module post_processing_mod
