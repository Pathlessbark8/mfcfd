module generate_connectivity_mod

	use data_structure_mod

	contains

	subroutine generate_connectivity()

		implicit none

		integer :: i, k
		real*8 :: nx, ny

		! do i = 1, local_points
		! 	nx = point(i)%nx
		! 	ny = point(i)%ny
		! 	if(point(i)%flag_1 == 1) then
        !         call get_interior_neighbours(i, nx, ny)
        !     else if(point(i)%flag_1 == 2) then
        !         call get_wall_boundary_neighbours(i, nx, ny)
        !     else if(point(k)%flag_1 == 3) then
		! 		call get_outer_boundary_neighbours(i, nx, ny)
        !     end if
		! enddo
		do k = 1, interior_points
					i = interior_points_index(k)
					nx = point(i)%nx
					ny = point(i)%ny
					call get_interior_neighbours(i, nx, ny)
		enddo

		do k = 1, wall_points
					i = wall_points_index(k)
					nx = point(i)%nx
					ny = point(i)%ny
					call get_wall_boundary_neighbours(i, nx, ny)
		enddo

		do k = 1, outer_points
					i = outer_points_index(k)
					nx = point(i)%nx
					ny = point(i)%ny
					call get_outer_boundary_neighbours(i, nx, ny)
		enddo

	end subroutine 

	subroutine get_interior_neighbours(i, nx, ny)

		implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
		integer :: i, r, count, nbh

			xi = point(i)%x
			yi = point(i)%y

			tx = ny
			ty = -nx

			point(i)%xpos_nbhs = 0
			point(i)%xneg_nbhs = 0
			point(i)%ypos_nbhs = 0
			point(i)%yneg_nbhs = 0

			do r=1, point(i)%nbhs
				nbh = point(i)%conn(r)
				xk = point(nbh)%x
				yk = point(nbh)%y

				delx = xk - xi
				dely = yk - yi
	!
				dels = delx*tx + dely*ty
				deln = delx*nx + dely*ny
	!
				if(dels .le. 0.0) then
							
					point(i)%xpos_nbhs = point(i)%xpos_nbhs + 1;

					count = point(i)%xpos_nbhs;
					point(i)%xpos_conn(count) = nbh;

				endif

				if(dels .gt. 0.0) then
					
					point(i)%xneg_nbhs = point(i)%xneg_nbhs + 1;

					count = point(i)%xneg_nbhs;
					point(i)%xneg_conn(count) = nbh;

				endif
							
				if(deln .le. 0.0) then
							
					point(i)%ypos_nbhs = point(i)%ypos_nbhs + 1;

					count = point(i)%ypos_nbhs;
					point(i)%ypos_conn(count) = nbh;

				endif
				
				if(deln .gt. 0.0) then
						
					point(i)%yneg_nbhs = point(i)%yneg_nbhs + 1;

					count = point(i)%yneg_nbhs;
					point(i)%yneg_conn(count) = nbh;

				endif
							
			enddo		

	end subroutine				

subroutine get_wall_boundary_neighbours(i, nx, ny)

		implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
		integer :: i, r, count, nbh
					

				xi = point(i)%x
				yi = point(i)%y

				tx = ny
                ty = -nx

				point(i)%xpos_nbhs = 0
                point(i)%xneg_nbhs = 0
                point(i)%yneg_nbhs = 0

				do r=1, point(i)%nbhs
				!
						nbh = point(i)%conn(r)

                        xk = point(nbh)%x
                        yk = point(nbh)%y
                        
                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0) then
						!
                                point(i)%xpos_nbhs = point(i)%xpos_nbhs + 1;

                                count = point(i)%xpos_nbhs;
                                point(i)%xpos_conn(count) = nbh;

						endif
						!
						if(dels .gt. 0.0) then
						!
                                point(i)%xneg_nbhs = point(i)%xneg_nbhs + 1;

                                count = point(i)%xneg_nbhs;
                                point(i)%xneg_conn(count) = nbh;
                                
						endif
						!
						point(i)%yneg_nbhs = point(i)%yneg_nbhs + 1;

						count = point(i)%yneg_nbhs;
						point(i)%yneg_conn(count) = nbh;

				enddo		

	end subroutine				
!
!
!
subroutine get_outer_boundary_neighbours(i, nx, ny)
!
!
		implicit none
!
		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
		integer :: i, r, count, nbh
!					
!
				xi = point(i)%x
				yi = point(i)%y

				tx = ny
                ty = -nx
                
                point(i)%xpos_nbhs = 0
                point(i)%xneg_nbhs = 0
                point(i)%ypos_nbhs = 0

				do r=1, point(i)%nbhs
				!
						nbh = point(i)%conn(r)

                        xk = point(nbh)%x
                        yk = point(nbh)%y
                        
                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0) then
						!
                                point(i)%xpos_nbhs = point(i)%xpos_nbhs + 1;

                                count = point(i)%xpos_nbhs;
                                point(i)%xpos_conn(count) = nbh;

						endif
						!
						if(dels .gt. 0.0) then
						!
                                point(i)%xneg_nbhs = point(i)%xneg_nbhs + 1;

                                count = point(i)%xneg_nbhs;
                                point(i)%xneg_conn(count) = nbh;

						endif
						!
						
						point(i)%ypos_nbhs = point(i)%ypos_nbhs + 1;

						count = point(i)%ypos_nbhs;
						point(i)%ypos_conn(count) = nbh;

						
				enddo		
!
!
	end subroutine				
!
!
!					
end module 		
						
						
						
						
						
			
