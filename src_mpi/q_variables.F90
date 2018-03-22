module q_variables_mod
!
!
	use data_structure_mod
!
contains
!
!
!
		subroutine eval_q_variables()
!
			implicit none
				
				integer :: k
				real*8 :: rho, u1, u2, pr, beta
				real*8 :: two_times_beta
!

 				do k = 1, local_points
				!
						rho = point(k)%rho
						u1 = point(k)%u1
						u2 = point(k)%u2
						pr = point(k)%pr			

		                beta = 0.5*rho/pr

						point(k)%q(1) = dlog(rho) + (dlog(beta)*2.5) - beta*(u1*u1 + u2*u2)

						two_times_beta = 2.0d0*beta

						point(k)%q(2) = two_times_beta*u1
						
						point(k)%q(3) = two_times_beta*u2
						
						point(k)%q(4) = -two_times_beta												
				!
				enddo
!				
		end subroutine 	
!
!
!
		subroutine eval_q_derivatives()
!
!
			implicit none
			!
			
				integer :: i, k, r, nbh
				
				real*8 :: x_i, y_i, x_k, y_k
				
				real*8 :: delx, dely, dist, weights
				real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
				real*8 :: sum_delx_delq(4), sum_dely_delq(4)
				real*8 :: det, delq, temp
				real*8 :: one_by_det
!				
!				
				do i = 1, local_points
				!
						
						x_i = point(i)%x
						y_i = point(i)%y
						
						sum_delx_sqr = 0.d0
						sum_dely_sqr = 0.d0
						sum_delx_dely = 0.d0
!						
						sum_delx_delq = 0.d0
						sum_dely_delq = 0.d0
!						
						do k = 1, point(i)%nbhs
						!
							nbh = point(i)%conn(k)
							
							x_k = point(nbh)%x
							y_k = point(nbh)%y
							
							delx = x_k - x_i
							dely = y_k - y_i
							
							dist = dsqrt(delx*delx + dely*dely)
							weights = dist**power
	
							sum_delx_sqr = sum_delx_sqr + delx*delx*weights
							sum_dely_sqr = sum_dely_sqr + dely*dely*weights
!
							sum_delx_dely = sum_delx_dely + delx*dely*weights
!							
							sum_delx_delq = sum_delx_delq + weights*delx*(point(nbh)%q - point(i)%q)
							sum_dely_delq = sum_dely_delq + weights*dely*(point(nbh)%q - point(i)%q)
!						!
						enddo							
!						
						det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
						one_by_det = 1.0d0/det
!									
						point(i)%qx = (sum_delx_delq*sum_dely_sqr - sum_dely_delq*sum_delx_dely)*one_by_det
						point(i)%qy = (sum_dely_delq*sum_delx_sqr - sum_delx_delq*sum_delx_dely)*one_by_det
!
!
!                                call test(point(i)%x,point(i)%y,point(i)%qx,0)
				enddo	
!
!
		end subroutine 	
!
!
!
		subroutine qtilde_to_primitive(qtilde, u1, u2, rho, pr)
!		!
			implicit none
!			!
				real*8 :: qtilde(4), u1, u2, rho, pr
				real*8 :: beta, temp, temp1, temp2
				real*8 :: q1, q2, q3, q4
!				
					q1 = qtilde(1)
					q2 = qtilde(2)
					q3 = qtilde(3)
					q4 = qtilde(4)
!				
					beta = -q4*0.5d0
!
					temp = 0.5d0/beta
!
			        u1 = q2*temp
			        u2 = q3*temp
!
			        temp1 = q1 + beta*(u1*u1 + u2*u2)
			        temp2 = temp1 - (dlog(beta)/(gamma-1))
!
					rho = exp(temp2)
					pr = rho*temp
!					
		end subroutine 		
!				

!
!
end module q_variables_mod
