module generate_connectivity_mod
        
        use data_structure_mod

        contains

        subroutine generate_connectivity()

                implicit none

                integer :: i, k
		real*8 :: nx, ny

                do k = 1, interior_points
                                        i = interior_points_index(k)
                                        nx = p%nx(i)
                                        ny = p%ny(i)
                                        call get_interior_neighbours(i, nx, ny)
                enddo

                do k = 1, wall_points
                                        i = wall_points_index(k)
                                        nx = p%nx(i)
                                        ny = p%ny(i)
                                        call get_wall_boundary_neighbours(i, nx, ny)
                enddo

                do k = 1, outer_points
                                        i = outer_points_index(k)
                                        nx = p%nx(i)
                                        ny = p%ny(i)
                                        call get_outer_boundary_neighbours(i, nx, ny)
                enddo

        end subroutine 

        subroutine get_interior_neighbours(i, nx, ny)

                implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
                integer :: i, r, count, nbh

                xi = p%x(i)
                yi = p%y(i)

                tx = ny
                ty = -nx

                p%xpos_nbhs(i) = 0
                p%xneg_nbhs(i) = 0
                p%ypos_nbhs(i) = 0
                p%yneg_nbhs(i) = 0

                do r=1, p%nbhs(i)
                        nbh = p%conn(i,r)
                        xk = p%x(nbh)
                        yk = p%y(nbh)

                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0) then

                                p%xpos_nbhs(i) = p%xpos_nbhs(i) + 1;

                                count = p%xpos_nbhs(i);
                                p%xpos_conn(i,count) = nbh;

                        endif

                        if(dels .gt. 0.0) then
                        
                                p%xneg_nbhs(i) = p%xneg_nbhs(i) + 1;

                                count = p%xneg_nbhs(i);
                                p%xneg_conn(i,count) = nbh;

                        endif

                        if(deln .le. 0.0) then

                                p%ypos_nbhs(i) = p%ypos_nbhs(i) + 1;

                                count = p%ypos_nbhs(i);
                                p%ypos_conn(i,count) = nbh;

                        endif
        
                        if(deln .gt. 0.0) then

                                p%yneg_nbhs(i) = p%yneg_nbhs(i) + 1;

                                count = p%yneg_nbhs(i);
                                p%yneg_conn(i,count) = nbh;
        
                        endif
                                

                enddo

        end subroutine

subroutine get_wall_boundary_neighbours(i, nx, ny)

                implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
                integer :: i, r, count, nbh


                xi = p%x(i)
                yi = p%y(i)

                tx = ny
                ty = -nx

                p%xpos_nbhs(i) = 0
                p%xneg_nbhs(i) = 0
                p%yneg_nbhs(i) = 0

                do r=1, p%nbhs(i)

                        nbh = p%conn(i,r)

                        xk = p%x(nbh)
                        yk = p%y(nbh)

                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0) then

                                p%xpos_nbhs(i) = p%xpos_nbhs(i) + 1;

                                count = p%xpos_nbhs(i);
                                p%xpos_conn(i,count) = nbh;

                        endif

                        if(dels .gt. 0.0) then

                                p%xneg_nbhs(i) = p%xneg_nbhs(i) + 1;

                                count = p%xneg_nbhs(i);
                                p%xneg_conn(i,count) = nbh;
  
                        endif

                        p%yneg_nbhs(i) = p%yneg_nbhs(i) + 1;

                        count = p%yneg_nbhs(i);
                        p%yneg_conn(i,count) = nbh;

                enddo

        end subroutine


subroutine get_outer_boundary_neighbours(i, nx, ny)

                implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
                integer :: i, r, count, nbh


                xi = p%x(i)
                yi = p%y(i)

                tx = ny
                ty = -nx
                
                p%xpos_nbhs(i) = 0
                p%xneg_nbhs(i) = 0
                p%ypos_nbhs(i) = 0

                do r=1, p%nbhs(i)

                        nbh = p%conn(i,r)

                        xk = p%x(nbh)
                        yk = p%y(nbh)
  
                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0) then

                                p%xpos_nbhs(i) = p%xpos_nbhs(i) + 1;

                                count = p%xpos_nbhs(i);
                                p%xpos_conn(i,count) = nbh;

                        endif

                        if(dels .gt. 0.0) then

                                p%xneg_nbhs(i) = p%xneg_nbhs(i) + 1;

                                count = p%xneg_nbhs(i);
                                p%xneg_conn(i,count) = nbh;

                        endif


                        p%ypos_nbhs(i) = p%ypos_nbhs(i) + 1;

                        count = p%ypos_nbhs(i);
                        p%ypos_conn(i,count) = nbh;


                enddo

        end subroutine

end module
