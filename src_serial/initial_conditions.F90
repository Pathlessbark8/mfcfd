module initial_conditions_mod

        use data_structure_mod
        use parameter_mod

contains

        subroutine initial_conditions()

                implicit none
                
                integer :: k,i,j, dummy, cycles
                ! integer :: nx=320, ny=120

                real*8 :: phi1_temp1, phi1_temp2, phi1_temp3, phi1_temp4
                real*8 :: phi2_temp1, phi2_temp2, phi2_temp3, phi2_temp4


                call setup_case_parameters()

                open(unit=109, file="vor_area", form="formatted", action="read")
                

                write(*,*) '%%%%%%%%%%%%-Reading Voronoi Area-%%%%%%%%%%%%'
                write(*,*) 

                do k=1, max_points
                        read(109,*) point%vor_area(k)
                end do

                close(unit=109)

                if(solution_restart .eq. 0) then

                        do k=1, max_points
                                point%prim(1,k) = q_init(1)
                                point%prim(2,k) = q_init(2)
                                point%prim(3,k) = q_init(3)
                                point%prim(4,k) = q_init(4)
                        enddo
                else

                        write(*,*) '%%%%%%%%%%%%-Using Restart file-%%%%%%%%%%%%'
                        write(*,*) 

                        open(unit=105, file="restart.dat", form='formatted', action="read")
                        read(105,*) dummy, itr, res_old
                        do k=1, max_points
                               read(105,*) dummy, dummy, dummy, dummy, &
                                       & point%prim(1,k), point%prim(2,k), point%prim(3,k), &
                                       & point%prim(4,k), dummy, dummy, dummy
                        end do

                        close(unit=105)

                endif

               if(phi_load .eq. 1) then

                        write(*,*) '%%%%%%%%%%%%-Loading Phi Vectors for Phi1 and Phi2 Base-%%%%%%%%%%%%'
                        write(*,*) 

                        open(unit=105, file="phi_vector.dat", form='formatted', action="read")
                        do k=1, max_points
                                read(105,*) point%phi1(1,k), point%phi1(2,k), point%phi1(3,k), &
                                & point%phi1(4,k), point%phi2(1,k), point%phi2(2,k), point%phi2(3,k), &
                                & point%phi2(4,k)
                        end do

                        close(unit=105)

                endif

                if(phi_load .eq. 2) then

                        write(*,*) '%%%%%%%%%%%%-Loading Phi Vectors for Phi1 and Phi2 with Epsilon-%%%%%%%%%%%%'
                        write(*,*) 
                        write(*,*) 'Phi Epsilon value is set to: ', phi_epsilon

                        open(unit=105, file="phi_vector.dat", form='formatted', action="read")
                        open(unit=106, file="phi_vector_d.dat", form='formatted', action="read")
                        open(unit=107, file="phi_vector_diff.dat", form='formatted', status="replace", action="write")

                        do k=1, max_points
                                read(105,*) point%phi1(1,k), point%phi1(2,k), point%phi1(3,k), &
                                & point%phi1(4,k), point%phi2(1,k), point%phi2(2,k), point%phi2(3,k), &
                                & point%phi2(4,k)
                                read(106,*) phi1_temp1, phi1_temp2, phi1_temp3, phi1_temp4, phi2_temp1, phi2_temp2, phi2_temp3, phi2_temp4
                               
                                point%phi1(1,k) = point%phi1(1,k) + (phi_epsilon * phi1_temp1)
                                point%phi1(2,k) = point%phi1(2,k) + (phi_epsilon * phi1_temp2)
                                point%phi1(3,k) = point%phi1(3,k) + (phi_epsilon * phi1_temp3)
                                point%phi1(4,k) = point%phi1(4,k) + (phi_epsilon * phi1_temp4)

                                point%phi2(1,k) = point%phi2(1,k) + (phi_epsilon * phi2_temp1)
                                point%phi2(2,k) = point%phi2(2,k) + (phi_epsilon * phi2_temp2)
                                point%phi2(3,k) = point%phi2(3,k) + (phi_epsilon * phi2_temp3)
                                point%phi2(4,k) = point%phi2(4,k) + (phi_epsilon * phi2_temp4)

                                do i=1, 4
                                        if(point%phi1(i,k) .gt. 1.0d0) then
                                                point%phi1(i,k) = 1.0d0
                                        endif
                                        if(point%phi2(i,k) .gt. 1.0d0) then
                                                point%phi2(i,k) = 1.0d0
                                        endif
                                        if(point%phi1(i,k) .lt. 0.0d0) then
                                                point%phi1(i,k) = 0.0d0
                                        endif
                                        if(point%phi2(i,k) .lt. 0.0d0) then
                                                point%phi2(i,k) = 0.0d0
                                        endif
                                enddo

                                write(107,'(8e30.20)') point%phi1(1,k), point%phi1(2,k), point%phi1(3,k), point%phi1(4,k), point%phi2(1,k), point%phi2(2,k), point%phi2(3,k), point%phi2(4,k)
                        end do

                        close(unit=105)
                        close(unit=106)
                        close(unit=107)

                endif

                open(unit=110, file="phi_vector_run.dat", action="write", status="replace")
                do i=1, max_points
                        write(110,'(8e30.20)') point%phi1(1,i), point%phi1(2,i), point%phi1(3,i), point%phi1(4,i), point%phi2(1,i), point%phi2(2,i), point%phi2(3,i), point%phi2(4,i)
                enddo
                close(unit=110)
        end subroutine

end module initial_conditions_mod
