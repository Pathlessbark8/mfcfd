module post_processing_mod

        use data_structure_mod

contains

        subroutine print_primal_output()


                implicit none

                integer :: i
                real*8 :: sos, vel_mag, mach_number
                character(len=64) :: sfile

                sfile = 'output.dat'

                OPEN(UNIT=501,FILE=trim(sfile))
                
                write(501,*)max_points, it, res_old

                do i = 1, max_points

                        sos = dsqrt(gamma*point%prim(4,i)/point%prim(1,i))
                        vel_mag = point%prim(2,i)*point%prim(2,i) + point%prim(3,i)*point%prim(3,i)

                        mach_number = dsqrt(vel_mag/sos)

                        write(501,'(2e30.20,2i4,7e30.20)')point%x(i), point%y(i), point%flag_1(i), point%qtdepth(i),&
                                & point%prim(1,i),point%prim(2,i),point%prim(3,i),&
                                & point%prim(4,i), mach_number, point%entropy(i), point%sensor(i)
                end do

                close(501)

        end subroutine

end module post_processing_mod