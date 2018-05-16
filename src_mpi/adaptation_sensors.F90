module adaptation_sensors_mod
!
	use data_structure_mod
!
	contains
!
!
		subroutine sensor_D2_distance()
!
			implicit none
!
			integer :: i, j, r
			real*8 :: u1_i, u2_i, pr_i, rho_i, x_i, y_i
			real*8 :: u1_j, u2_j, pr_j, rho_j, x_j, y_j
			
			real*8 :: max, min, dist_ij, u_sqr
			real*8 :: temp1, temp2, temp3
			real*8 :: D2_dist, D2
						
!
			do i = 1, local_points
			!
				u1_i = point%prim(2,i)
				u2_i = point%prim(3,i)
				rho_i = point%prim(1,i)
				pr_i = point%prim(4,i)
				x_i = point%x(i)
				y_i = point%y(i)
			!
				max = 0.0
				min = 1.0e10
				
				do r = 1, point%nbhs(i)
					j = point%conn(i,r)
!
					u1_j = point%prim(2,i)
					u2_j = point%prim(3,i)
					pr_j = point%prim(1,i)
					rho_j = point%prim(4,i)
					x_j = point%x(i)
					y_j = point%y(i)
					
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
					
					D2 = dabs(D2_dist/dist_ij)
					
					if(D2 > max) then 
						max = D2
					endif
					
					if(D2 < min) then
						min = D2
					endif
					
!					print*, i, j, D2_dist, D2
				!
				enddo
!
				point%sensor(i) = max/min
				
			enddo		
!
!
		end subroutine sensor_D2_distance
!
!									
!
end module adaptation_sensors_mod
!		
!
!
