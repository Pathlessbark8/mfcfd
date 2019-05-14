module q_variables_mod

       
        use data_structure_mod

contains

                subroutine eval_q_variables()

                        implicit none

                                integer :: k
				real*8 :: rho, u1, u2, pr, beta
				real*8 :: two_times_beta
!

                                do k = 1, max_points

                                                rho = point%prim(1,k)
                                                u1 = point%prim(2,k)
                                                u2 = point%prim(3,k)
                                                pr = point%prim(4,k)

                                beta = 0.5d0*rho/pr

                                                point%q(1,k) = dlog(rho) + (dlog(beta)*2.5d0) - beta*(u1*u1 + u2*u2)

                                                two_times_beta = 2.0d0*beta

                                                point%q(2,k) = two_times_beta*u1

                                                point%q(3,k) = two_times_beta*u2

                                                point%q(4,k) = -two_times_beta
                                enddo

                        end subroutine 



                subroutine eval_q_derivatives()


                        implicit none


                                integer :: i, k, r, nbh

				real*8 :: x_i, y_i, x_k, y_k

				real*8 :: delx, dely, dist, weights
				real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
				real*8 :: sum_delx_delq(4), sum_dely_delq(4)
				real*8 :: det, delq, temp
				real*8 :: one_by_det

                                do i = 1, local_points


                                        x_i = point%x(i)
                                        y_i = point%y(i)

                                        sum_delx_sqr = 0.d0
                                        sum_dely_sqr = 0.d0
                                        sum_delx_dely = 0.d0

                                        sum_delx_delq = 0.d0
                                        sum_dely_delq = 0.d0

                                        point%qm(1, :, i) = point%q(:, i)
                                        point%qm(2, :, i) = point%q(:, i)

                                        do k = 1, point%nbhs(i)

                                                nbh = point%conn(i,k)

                                                do r = 1, 4
                                                        if(point%q(r,nbh) > point%qm(1,r,i)) then
                                                                point%qm(1,r,i) = point%q(r,nbh)
                                                        end if
                                                        if(point%q(r,nbh) < point%qm(2,r,i)) then
                                                                point%qm(2,r,i) = point%q(r,nbh)
                                                        end if
                                                end do

                                                x_k = point%x(nbh)
                                                y_k = point%y(nbh)

                                                delx = x_k - x_i
                                                dely = y_k - y_i

                                                dist = dsqrt(delx*delx + dely*dely)
                                                weights = dist**power

                                                sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                                                sum_dely_sqr = sum_dely_sqr + dely*dely*weights

                                                sum_delx_dely = sum_delx_dely + delx*dely*weights

                                                sum_delx_delq = sum_delx_delq + weights*delx*(point%q(:,nbh) - point%q(:,i))
                                                sum_dely_delq = sum_dely_delq + weights*dely*(point%q(:,nbh) - point%q(:,i))


                                        enddo

                                        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                                        one_by_det = 1.0d0/det

                                        point%dq(1,:,i) = (sum_delx_delq*sum_dely_sqr&
                                               & - sum_dely_delq*sum_delx_dely)*one_by_det
                                        point%dq(2,:,i) = (sum_dely_delq*sum_delx_sqr&
                                                       &- sum_delx_delq*sum_delx_dely)*one_by_det




                                enddo


                end subroutine 



                subroutine qtilde_to_primitive(qtilde, u1, u2, rho, pr)

                        implicit none

				real*8 :: qtilde(4), u1, u2, rho, pr
				real*8 :: beta, temp, temp1, temp2
				real*8 :: q1, q2, q3, q4

                                q1 = qtilde(1)
                                q2 = qtilde(2)
                                q3 = qtilde(3)
                                q4 = qtilde(4)
        
                                beta = -q4*0.5d0

                                temp = 0.5d0/beta

                                u1 = q2*temp
                                u2 = q3*temp

                                temp1 = q1 + beta*(u1*u1 + u2*u2)
                                temp2 = temp1 - (dlog(beta)/(gamma-1))

                                rho = dexp(temp2)
                                pr = rho*temp

                end subroutine 



end module q_variables_mod
