module point_normals_mod
!
	use data_structure_mod
!
	contains
	!
	!
	subroutine compute_normals()
!
                implicit none
!
                double precision :: lx, ly, mx, my, rx, ry
                double precision :: nx1, nx2, ny1, ny2, nx, ny
                double precision :: det
                !
                integer:: i, j, k, l, m, r
!
!	
!	Finding the normals for the points on the shapes ..   
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
						nx1 = my - ly
						nx2 = ry - my
!
						ny1 = mx - lx
						ny2 = rx - mx
!	
						nx = 0.5d0*(nx1 + nx2)
						ny = 0.5d0*(ny1 + ny2)
!
						det = dsqrt(nx*nx + ny*ny)
!
						nx = -nx/det
						ny = ny/det
!
						k = shape_points_index(j, i)
!
						point(k)%nx = nx
						point(k)%ny = ny
!
					enddo																		
				enddo		
!			
!	Finding the normals for the outer boundary points ..
!
				do i = 1, outer_points
!
						if(i .eq. 1) then 
!						!
								m = outer_points_index(i)
								r = outer_points_index(i+1)
								l = outer_points_index(1)
						endif
						!
						if(i .gt. 1 .and. i .lt. outer_points) then 
						!
								m = outer_points_index(i)
								r = outer_points_index(i+1)
								l = outer_points_index(i-1)
						endif
						!
						if(i .eq. outer_points) then 
!					!
								m = outer_points_index(i)
								r = outer_points_index(1)
								l = outer_points_index(i-1)
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
						nx1 = my - ly
						nx2 = ry - my
!
						ny1 = mx - lx
						ny2 = rx - mx
!	
						nx = 0.5d0*(nx1 + nx2)
						ny = 0.5d0*(ny1 + ny2)
!
						det = dsqrt(nx*nx + ny*ny)
!
						nx = -nx/det
						ny = ny/det
!
						k = outer_points_index(i)
!
						point(k)%nx = nx
						point(k)%ny = ny
!
					enddo		
!
!
!	The below lines of code are temporary. In future
!	when we want the normals of the interior points
!	to be taken into account then we need to evaluate 
!	the normals here ..
!
					do i = 1, interior_points

						k = interior_points_index(i)
!
						point(k)%nx = 0.0d0
						point(k)%ny = 1.0d0
!
					enddo						
!
!
!	Suppose the normals of the interior points are available
!	and we still want to do upwinding along the cartesian coordinates
!	then the following portion of the code ensures it ..
!
!
					if(interior_points_normal_flag .eq. 0) then
						do i = 1, interior_points
							k = interior_points_index(i)
							point(k)%nx = 0.d0
							point(k)%ny = 1.d0
						enddo
					endif									
!
!
        end subroutine 
!
end module 		
!
!						
