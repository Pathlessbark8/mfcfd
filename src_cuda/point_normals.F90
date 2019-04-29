module point_normals_mod

        use data_structure_mod

        contains


        subroutine compute_normals()

                implicit none

                double precision :: lx, ly, mx, my, rx, ry
                double precision :: nx1, nx2, ny1, ny2, nx, ny
                double precision :: det
                
                integer:: i, j, k, l, m, r


!Finding the normals for the points on the shapes ..   

                do i = 1, shape_points
                        m = shape_points_index(i)
                        l = point%left(m)
                        r = point%right(m)

                        
                        
                        lx = point%x(l)
                        ly = point%y(l)

                        mx = point%x(m)
                        my = point%y(m)

                        rx = point%x(r)
                        ry = point%y(r)

                        nx1 = my - ly
                        nx2 = ry - my

                        ny1 = mx - lx
                        ny2 = rx - mx

                        nx = 0.5d0*(nx1 + nx2)
                        ny = 0.5d0*(ny1 + ny2)

                        det = dsqrt(nx*nx + ny*ny)

                        nx = -nx/det
                        ny = ny/det

                        point%nx(m) = nx
                        point%ny(m) = ny

                enddo
                                

!	Finding the normals for the outer boundary points ..

                do i = 1, outer_points
                        m = outer_points_index(i)
                        l = point%left(m)
                        r = point%right(m)

                        lx = point%x(l)
                        ly = point%y(l)

                        mx = point%x(m)
                        my = point%y(m)

                        rx = point%x(r)
                        ry = point%y(r)
   
                        nx1 = my - ly
                        nx2 = ry - my

                        ny1 = mx - lx
                        ny2 = rx - mx

                        nx = 0.5d0*(nx1 + nx2)
                        ny = 0.5d0*(ny1 + ny2)

                        det = dsqrt(nx*nx + ny*ny)

                        nx = -nx/det
                        ny = ny/det

                        point%nx(m) = nx
                        point%ny(m) = ny

                enddo

        


!	The below lines of code are temporary. In future
!	when we want the normals of the interior points
!	to be taken into account then we need to evaluate 
!	the normals here ..
!
                do i = 1, interior_points

                        k = interior_points_index(i)

                        point%nx(k) = 0.0d0
                        point%ny(k) = 1.0d0

                enddo


!	Suppose the normals of the interior points are available
!	and we still want to do upwinding along the cartesian coordinates
!	then the following portion of the code ensures it ..


                if(interior_points_normal_flag .eq. 0) then
                        do i = 1, interior_points
                                k = interior_points_index(i)
                                point%nx(k) = 0.d0
                                point%ny(k) = 1.d0
                        enddo
                endif

        end subroutine 

end module 


