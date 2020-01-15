module generate_connectivity_mod
        
        use DATA_STRUCTURE_MOD_DIFF

        contains

        subroutine generate_connectivity()

                implicit none

                integer :: i, k
		real*8 :: nx, ny

                do k = 1, interior_points
                                        i = interior_points_index(k)
                                        nx = point%nx(i)
                                        ny = point%ny(i)
                                        call get_interior_neighbours(i, nx, ny)
                enddo

                do k = 1, wall_points
                                        i = wall_points_index(k)
                                        nx = point%nx(i)
                                        ny = point%ny(i)
                                        call get_wall_boundary_neighbours(i, nx, ny)
                enddo

                do k = 1, outer_points
                                        i = outer_points_index(k)
                                        nx = point%nx(i)
                                        ny = point%ny(i)
                                        call get_outer_boundary_neighbours(i, nx, ny)
                enddo

        end subroutine 

        subroutine get_interior_neighbours(i, nx, ny)

                implicit none

		real*8 :: xi, yi, xk, yk
		real*8 :: delx, dely, dels, deln
		real*8 :: nx, ny, tx, ty
                integer :: i, r, count, nbh

                xi = point%x(i)
                yi = point%y(i)

                tx = ny
                ty = -nx

                point%xpos_nbhs(i) = 0
                point%xneg_nbhs(i) = 0
                point%ypos_nbhs(i) = 0
                point%yneg_nbhs(i) = 0

                do r=1, point%nbhs(i)
                        nbh = point%conn(i,r)
                        xk = point%x(nbh)
                        yk = point%y(nbh)
                        
                        nbh = find_loc_f90(point%conn, 20, i, nbh)

                        delx = xk - xi
                        dely = yk - yi

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        if(dels .le. 0.0d0) then

                                point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;

                                count = point%xpos_nbhs(i);
                                point%xpos_conn(i,count) = nbh;

                        endif

                        if(dels .ge. 0.0d0) then
                        
                                point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;

                                count = point%xneg_nbhs(i);
                                point%xneg_conn(i,count) = nbh;

                        endif

                        if(deln .le. 0.0d0) then

                                point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1;

                                count = point%ypos_nbhs(i);
                                point%ypos_conn(i,count) = nbh;

                        endif
        
                        if(deln .ge. 0.0d0) then

                                point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1;

                                count = point%yneg_nbhs(i);
                                point%yneg_conn(i,count) = nbh;
        
                        endif

                enddo

                if(point%xpos_nbhs(i) == 0) then
                        print*,"xpos zero for interior point number:", i
                        stop
                elseif(point%xneg_nbhs(i) == 0) then
                        print*,"xneg zero for interior point number:", i
                        stop
                elseif(point%ypos_nbhs(i) == 0) then
                        print*,"ypos zero for interior point number:", i
                        stop
                elseif(point%yneg_nbhs(i) == 0) then
                        print*,"yneg zero for interior point number:", i
                        stop
                end if
                
        end subroutine

        subroutine get_wall_boundary_neighbours(i, nx, ny)
        
                        implicit none
        
        		real*8 :: xi, yi, xk, yk
        		real*8 :: delx, dely, dels, deln
        		real*8 :: nx, ny, tx, ty
                        integer :: i, r, count, nbh
        
        
                        xi = point%x(i)
                        yi = point%y(i)
        
                        tx = ny
                        ty = -nx
        
                        point%xpos_nbhs(i) = 0
                        point%xneg_nbhs(i) = 0
                        point%yneg_nbhs(i) = 0
        
                        do r=1, point%nbhs(i)
        
                                nbh = point%conn(i,r)
        
                                xk = point%x(nbh)
                                yk = point%y(nbh)

                                nbh = find_loc_f90(point%conn, 20, i, nbh)


                                delx = xk - xi
                                dely = yk - yi
        
                                dels = delx*tx + dely*ty
                                deln = delx*nx + dely*ny
        
                                if(dels .le. 0.0d0) then
        
                                        point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;
        
                                        count = point%xpos_nbhs(i);
                                        point%xpos_conn(i,count) = nbh;
        
                                endif
        
                                if(dels .ge. 0.0d0) then
        
                                        point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;
        
                                        count = point%xneg_nbhs(i);
                                        point%xneg_conn(i,count) = nbh;
          
                                endif
        
                                point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1;
        
                                count = point%yneg_nbhs(i);
                                point%yneg_conn(i,count) = nbh;
        
                        enddo
                        
                        if(point%xpos_nbhs(i) == 0) then
                                print*,"xpos zero for wall point number:", i
                                stop
                        elseif(point%xneg_nbhs(i) == 0) then
                                print*,"xneg zero for wall point number:", i
                                stop
                        elseif(point%yneg_nbhs(i) == 0) then
                                print*,"yneg zero for wall point number:", i
                                stop
                        end if
        
        end subroutine
        
        
        subroutine get_outer_boundary_neighbours(i, nx, ny)
        
                        implicit none
        
        		real*8 :: xi, yi, xk, yk
        		real*8 :: delx, dely, dels, deln
        		real*8 :: nx, ny, tx, ty
                        integer :: i, r, count, nbh
        
        
                        xi = point%x(i)
                        yi = point%y(i)
        
                        tx = ny
                        ty = -nx
                        
                        point%xpos_nbhs(i) = 0
                        point%xneg_nbhs(i) = 0
                        point%ypos_nbhs(i) = 0
        
                        do r=1, point%nbhs(i)
        
                                nbh = point%conn(i,r)
        
                                xk = point%x(nbh)
                                yk = point%y(nbh)

                                nbh = find_loc_f90(point%conn, 20, i, nbh)

          
                                delx = xk - xi
                                dely = yk - yi
        
                                dels = delx*tx + dely*ty
                                deln = delx*nx + dely*ny
        
                                if(dels .le. 0.0d0) then
        
                                        point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;
        
                                        count = point%xpos_nbhs(i);
                                        point%xpos_conn(i,count) = nbh;
        
                                endif
        
                                if(dels .ge. 0.0d0) then
        
                                        point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;
        
                                        count = point%xneg_nbhs(i);
                                        point%xneg_conn(i,count) = nbh;
        
                                endif
        
        
                                point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1;
        
                                count = point%ypos_nbhs(i);
                                point%ypos_conn(i,count) = nbh;
        
                        enddo
        
                        if(point%xpos_nbhs(i) == 0) then
                                print*,"xpos zero for outer point number:", i
                                stop
                        elseif(point%xneg_nbhs(i) == 0) then
                                print*,"xneg zero for outer point number:", i
                                stop
                        elseif(point%ypos_nbhs(i) == 0) then
                                print*,"ypos zero for outer point number:", i
                                stop
                        end if

        end subroutine

        integer function find_loc_f90(array, pointcount, pidx, nbhvalue)
                integer, dimension(:,:) :: array
                integer :: pointcount, pidx, nbhvalue, i
                do i = 1, pointcount
                        if (array(pidx, i) == nbhvalue) then
                                find_loc_f90 = i
                                return
                        endif
                end do
                print*, "warning could not find point in conn", pidx, nbhvalue
                stop
        end function

end module
