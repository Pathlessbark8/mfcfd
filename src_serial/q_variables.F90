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

                                do i = 1, max_points


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

                subroutine eval_q_double_derivatives()

                        implicit none

                        ! local variables
                        integer :: i
                        integer :: k, r, nbh
                        real*8 :: x_i, y_i, x_k, y_k
                        real*8 :: delx, dely, dist, weights
                        real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                        real*8 :: sum_delx_del_qx(4), sum_delx_del_qy(4), sum_dely_del_qx(4), sum_dely_del_qy(4)
                        real*8 :: det, delq, temp
                        real*8 :: one_by_det

                        do i=1, max_points

                                x_i = point%x(i)
                                y_i = point%y(i)

                                sum_delx_sqr = 0.d0
                                sum_dely_sqr = 0.d0
                                sum_delx_dely = 0.d0

                                sum_delx_del_qx = 0.d0
                                sum_delx_del_qy = 0.d0
                                sum_dely_del_qx = 0.d0
                                sum_dely_del_qy = 0.d0

                                do k = 1, point%nbhs(i)

                                        nbh = point%conn(i,k)
                                        
                                        x_k = point%x(nbh)
                                        y_k = point%y(nbh)

                                        delx = x_k - x_i
                                        dely = y_k - y_i

                                        dist = dsqrt(delx*delx + dely*dely)
                                        weights = dist**power

                                        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                                        sum_dely_sqr = sum_dely_sqr + dely*dely*weights

                                        sum_delx_dely = sum_delx_dely + delx*dely*weights

                                        sum_delx_del_qx = sum_delx_del_qx + weights*delx*(point%dq(1,:,nbh) - point%dq(1,:,i))
                                        sum_delx_del_qy = sum_delx_del_qy + weights*delx*(point%dq(2,:,nbh) - point%dq(2,:,i))
                                        sum_dely_del_qx = sum_dely_del_qx + weights*dely*(point%dq(1,:,nbh) - point%dq(1,:,i))
                                        sum_dely_del_qy = sum_dely_del_qy + weights*dely*(point%dq(2,:,nbh) - point%dq(2,:,i))

                                enddo

                                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                                one_by_det = 1.0d0/det

                                point%ddq(1,:,i) = (sum_delx_del_qx*sum_dely_sqr&
                                & - sum_dely_del_qx*sum_delx_dely)*one_by_det
                                point%ddq(2,:,i) = (sum_delx_del_qy*sum_dely_sqr&
                                & - sum_dely_del_qy*sum_delx_dely)*one_by_det
                                point%ddq(3,:,i) = (sum_dely_del_qy*sum_delx_sqr&
                                &- sum_delx_del_qy*sum_delx_dely)*one_by_det
                        end do

                end subroutine

                subroutine eval_q_inner_loop()
                        implicit none

                        integer :: i
                        integer :: k, r, nbh
                        real*8 :: x_i, y_i, x_k, y_k
                        real*8 :: delx, dely, dist, weights
                        real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                        real*8 :: det, temp
                        real*8 :: one_by_det
                        real*8 :: sum_delx_delq1, sum_delx_delq2, sum_delx_delq3, sum_delx_delq4
                        real*8 :: sum_dely_delq1, sum_dely_delq2, sum_dely_delq3, sum_dely_delq4
                        real*8 :: q1, q2, q3, q4
                
                        real*8 :: temp1, temp2

                        do i = 1, max_points
        
                                x_i = point%x(i)
                                y_i = point%y(i)
                
                                sum_delx_sqr = 0.d0
                                sum_dely_sqr = 0.d0
                                sum_delx_dely = 0.d0
                
                                temp1 = 0.d0
                                temp2 = 0.d0
                
                                sum_delx_delq1 = 0.d0
                                sum_delx_delq2 = 0.d0
                                sum_delx_delq3 = 0.d0
                                sum_delx_delq4 = 0.d0
                
                                sum_dely_delq1 = 0.d0
                                sum_dely_delq2 = 0.d0
                                sum_dely_delq3 = 0.d0
                                sum_dely_delq4 = 0.d0
                
                                q1 = point%q(1, i)
                                q2 = point%q(2, i)
                                q3 = point%q(3, i)
                                q4 = point%q(4, i)
                
                                do k = 1, point%nbhs(i)
                
                                        nbh = point%conn(i,k)
                                        
                                        x_k = point%x(nbh)
                                        y_k = point%y(nbh)
                
                                        delx = x_k - x_i
                                        dely = y_k - y_i
                
                                        dist = dsqrt(delx*delx + dely*dely)
                                        weights = dist**power
                
                                        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                                        sum_dely_sqr = sum_dely_sqr + dely*dely*weights
                
                                        sum_delx_dely = sum_delx_dely + delx*dely*weights
                
                                        temp1 = q1 - 0.5d0 * (delx * point%dq(1,1,i) + dely * point%dq(2,1,i)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,1,i) + 2.0d0*delx*dely * point%ddq(2,1,i) * dely*dely * point%ddq(3,1,i))
                                        temp2 = point%q(1,nbh) - 0.5d0 * (delx * point%dq(1,1,nbh) + dely * point%dq(2,1,nbh)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,1,nbh) + 2*delx*dely * point%ddq(2,1,nbh) * dely*dely * point%ddq(3,1,nbh))
                                        sum_delx_delq1 = sum_delx_delq1 + (weights * delx * (temp2 - temp1))
                                        sum_dely_delq1 = sum_dely_delq1 + (weights * dely * (temp2 - temp1))
                
                                        temp1 = q2 - 0.5d0 * (delx * point%dq(1,2,i) + dely * point%dq(2,2,i)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,2,i) + 2.0d0*delx*dely * point%ddq(2,2,i) * dely*dely * point%ddq(3,2,i))
                                        temp2 = point%q(2,nbh) - 0.5d0 * (delx * point%dq(1,2,nbh) + dely * point%dq(2,2,nbh)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,2,nbh) + 2*delx*dely * point%ddq(2,2,nbh) * dely*dely * point%ddq(3,2,nbh))
                                        sum_delx_delq2 = sum_delx_delq2 + (weights * delx * (temp2 - temp1))
                                        sum_dely_delq2 = sum_dely_delq2 + (weights * dely * (temp2 - temp1))
                
                                        temp1 = q3 - 0.5d0 * (delx * point%dq(1,3,i) + dely * point%dq(2,3,i)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,3,i) + 2.0d0*delx*dely * point%ddq(2,3,i) * dely*dely * point%ddq(3,3,i))
                                        temp2 = point%q(3,nbh) - 0.5d0 * (delx * point%dq(1,3,nbh) + dely * point%dq(2,3,nbh)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,3,nbh) + 2*delx*dely * point%ddq(2,3,nbh) * dely*dely * point%ddq(3,3,nbh))
                                        sum_delx_delq3 = sum_delx_delq3 + (weights * delx * (temp2 - temp1))
                                        sum_dely_delq3 = sum_dely_delq3 + (weights * dely * (temp2 - temp1))
                
                                        temp1 = q4 - 0.5d0 * (delx * point%dq(1,4,i) + dely * point%dq(2,4,i)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,4,i) + 2.0d0*delx*dely * point%ddq(2,4,i) * dely*dely * point%ddq(3,4,i))
                                        temp2 = point%q(4,nbh) - 0.5d0 * (delx * point%dq(1,4,nbh) + dely * point%dq(2,4,nbh)) + (1.0d0/12.0d0) * (delx*delx * point%ddq(1,4,nbh) + 2*delx*dely * point%ddq(2,4,nbh) * dely*dely * point%ddq(3,4,nbh))
                                        sum_delx_delq4 = sum_delx_delq4 + (weights * delx * (temp2 - temp1))
                                        sum_dely_delq4 = sum_dely_delq4 + (weights * dely * (temp2 - temp1))
                
                                enddo
                
                                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                                one_by_det = 1.0d0/det
                
                                point%temp(1,1,i) = one_by_det * (sum_delx_delq1 * sum_dely_sqr - sum_dely_delq1 * sum_delx_dely)
                                point%temp(1,2,i) = one_by_det * (sum_delx_delq2 * sum_dely_sqr - sum_dely_delq2 * sum_delx_dely)
                                point%temp(1,3,i) = one_by_det * (sum_delx_delq3 * sum_dely_sqr - sum_dely_delq3 * sum_delx_dely)
                                point%temp(1,4,i) = one_by_det * (sum_delx_delq4 * sum_dely_sqr - sum_dely_delq4 * sum_delx_dely)
                                point%temp(2,1,i) = one_by_det * (sum_dely_delq1 * sum_delx_sqr - sum_delx_delq1 * sum_delx_dely)
                                point%temp(2,2,i) = one_by_det * (sum_dely_delq2 * sum_delx_sqr - sum_delx_delq2 * sum_delx_dely)
                                point%temp(2,3,i) = one_by_det * (sum_dely_delq3 * sum_delx_sqr - sum_delx_delq3 * sum_delx_dely)
                                point%temp(2,4,i) = one_by_det * (sum_dely_delq4 * sum_delx_sqr - sum_delx_delq4 * sum_delx_dely)
        
                        enddo

                end subroutine

                subroutine eval_update_innerloop()
                        do i=1,max_points
                                point%dq(1,1,i) = point%temp(1,1,i)
                                point%dq(1,2,i) = point%temp(1,2,i)
                                point%dq(1,3,i) = point%temp(1,3,i)
                                point%dq(1,4,i) = point%temp(1,4,i)
                                point%dq(2,1,i) = point%temp(2,1,i)
                                point%dq(2,2,i) = point%temp(2,2,i)
                                point%dq(2,3,i) = point%temp(2,3,i)
                                point%dq(2,4,i) = point%temp(2,4,i)
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
