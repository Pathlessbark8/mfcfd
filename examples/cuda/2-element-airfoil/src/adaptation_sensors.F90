module adaptation_sensors_mod

        use data_structure_mod

        contains

        subroutine compute_adapt_sensor()


                implicit none

                integer :: i
                real*8 :: max_sensor

                max_sensor = 0.0d0                                

                call sensor_D2_distance()

                do i = 1, max_points
                        point%sensor(i) = point%D2_dist(i)
                enddo

                do i = 1, max_points
                        if(point%sensor(i) .ge. max_sensor) then
                                max_sensor = point%sensor(i)
                        endif
                enddo

                do i = 1, max_points
                        point%sensor(i) = point%sensor(i)/max_sensor
                enddo


        end subroutine 



        subroutine sensor_D2_distance()

                implicit none

                integer :: i, j, r
		real*8 :: u1_i, u2_i, pr_i, rho_i, x_i, y_i
		real*8 :: u1_j, u2_j, pr_j, rho_j, x_j, y_j

		real*8 :: maxi, mini, dist_ij, u_sqr
		real*8 :: temp1, temp2, temp3
		real*8 :: D2_dist1, D2_dist2, D2_dist
		real*8 :: max_D2
                max_D2 = 0.0d0

                do i = 1, max_points

                        u1_i = point%prim(2,i)
                        u2_i = point%prim(3,i)
                        rho_i = point%prim(1,i)
                        pr_i = point%prim(4,i)
                        x_i = point%x(i)
                        y_i = point%y(i)

                        maxi = 0.0
                        mini = 1.0e10

                        do r = 1, point%nbhs(i)
                                j = point%conn(i,r)

                                u1_j = point%prim(2,j)
                                u2_j = point%prim(3,j)
                                rho_j = point%prim(1,j)
                                pr_j = point%prim(4,j)
                                x_j = point%x(j)
                                y_j = point%y(j)

                                dist_ij = (x_j - x_i)**2 + (y_j - y_i)**2

                                temp1 = (pr_i*pr_i*rho_j)/(pr_j*pr_j*rho_i)
                                temp1 = dlog(temp1)
                                temp1 = (rho_j - rho_i)*temp1

                                u_sqr = (u1_j - u1_i)**2 + (u2_j - u2_i)**2
                                temp2 = (rho_j/rho_i) + (pr_i*rho_j)/(pr_j*rho_i)
                                temp2 = temp2*u_sqr
                                temp2 = (0.5*rho_i*rho_i/pr_i)*temp2

                                temp3 = rho_i*(((pr_i*rho_j)/(pr_j*rho_i)) - 1)
                                temp3 = temp3 + rho_j*(((pr_j*rho_i)/(pr_i*rho_j)) - 1) 
                                temp3 = 2*temp3

                                D2_dist = temp1 + temp2 + temp3 

                                if(D2_dist > maxi) then 
                                        maxi = D2_dist
                                endif

                        enddo

                        point%D2_dist(i) = maxi

                enddo

        end subroutine sensor_D2_distance

end module adaptation_sensors_mod
