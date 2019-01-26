module post_processing_mod


        use data_structure_mod


contains


        subroutine print_primal_output()


                implicit none

                integer :: i
                character(len=64) :: sfile
                character(len=10) :: itos

                sfile = 'sol.dat'

                OPEN(UNIT=501,FILE=trim(sfile))
                
                write(501,*)max_points

                do i = 1, max_points
                        write(501,*)point%global_id(i),point%flag_1(i),point%flag_2(i),point%x(i),&
                                point%y(i),point%prim(1,i),point%prim(2,i),point%prim(3,i),&
                                point%prim(4,i)
                end do

                close(501)



        end subroutine

end module post_processing_mod
