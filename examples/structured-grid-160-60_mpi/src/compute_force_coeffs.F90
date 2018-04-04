module compute_force_coeffs_mod
!
	use data_structure_mod
!
	contains
!
!
!
		subroutine	compute_cl_cd_cm()
!
!
			implicit none
!
			integer :: i, j, k
			integer :: l, m, r
			real*8 :: cp, temp
			real*8 :: lx, ly, mx, my, rx, ry
			real*8 :: ds1, ds2, ds
!
			real*8 :: H, V, pitch_mom
			real*8 :: nx, ny
!			
!
			OPEN(UNIT=201,FILE="cp-file.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!
!			
			temp = 0.5*rho_inf*Mach*Mach
!			
			H = 0.d0
			V = 0.d0
			pitch_mom = 0.0d0
!
			do j = 1, shapes
				do i = 1, shape_points(j)
					if(i .eq. 1) then 
					!
							m = shape_points_index(j, i) 				
							r = shape_points_index(j, i+1) 
							l = shape_points_index(j, 1) 	
					endif
					!
					if(i .gt. 1 .and. i .lt. shape_points(j)) then 
					!
							m = shape_points_index(j, i) 				
							r = shape_points_index(j, i+1) 				
							l = shape_points_index(j, i-1) 				
					endif
					!
					if(i .eq. shape_points(j)) then 
					!
							m = shape_points_index(j, i) 				
							r = shape_points_index(j, 1)
							l = shape_points_index(j, i-1) 
					endif
!
					lx = point(l)%x
					ly = point(l)%y
!							
					mx = point(m)%x
					my = point(m)%y
!							
					rx = point(r)%x
					ry = point(r)%y
!
					ds1 = (mx - lx)**2 + (my - ly)**2
					ds1 = dsqrt(ds1)
!							
					ds2 = (rx - mx)**2 + (ry - my)**2
					ds2 = dsqrt(ds2)
!							
					ds = 0.5*(ds1 + ds2)				
!
					k = shape_points_index(j, i)
					nx = point(k)%nx
					ny = point(k)%ny
!					
					cp = point(k)%pr - pr_inf
					cp = -cp/temp
!
					write(201, *) point(k)%x, cp
!				
					H = H + cp*nx*ds
					V = V + cp*ny*ds					
					pitch_mom = pitch_mom + (-cp*ny*ds*(mx - 0.25) + cp*nx*ds*(my))
!					
				enddo
			enddo
!
			Cl = V*dcos(theta) - H*dsin(theta)
			Cd = H*dcos(theta) + V*dsin(theta)
			Cm = pitch_mom
!
			CLOSE(UNIT=201)	
!
!
		end subroutine 
!
end module compute_force_coeffs_mod
!		
!
!
