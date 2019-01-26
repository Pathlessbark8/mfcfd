module point_preprocessor_mod

        use data_structure_mod

contains

        subroutine read_input_point_data()

                implicit none

                integer:: i, k, r
                integer :: wall_temp,outer_temp,interior_temp
                character(len=64) :: part_grid

                print*,'Reading points'
                part_grid = 'partGrid'

                OPEN(UNIT=101,FILE=trim(part_grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

                !TODO: Add asserts

                read(101,*) max_points
                allocate(point%x(max_points))
                allocate(point%y(max_points))
                allocate(point%local_id(max_points))
                allocate(point%global_id(max_points))
                allocate(point%flag_1(max_points))
                allocate(point%flag_2(max_points))
                allocate(point%nbhs(max_points))
                allocate(point%conn(max_points,15))
                allocate(point%nx(max_points))
                allocate(point%ny(max_points))
                allocate(point%left(max_points))
                allocate(point%right(max_points))
                
                wall_points = 0
                interior_points = 0
                outer_points = 0

                if(old_format == 0) then
                        do k = 1, max_points
        
                                read(101,*) point%local_id(k),point%global_id(k),point%x(k),&
                                & point%y(k),point%left(k),point%right(k), point%flag_1(k),point%flag_2(k),point%nbhs(k),&
                                & (point%conn(k,r),r=1,point%nbhs(k))
                                
                        !Storing the count for the point types
                                if(point%flag_1(k) == 0) then
                                        wall_points = wall_points + 1
                                else if(point%flag_1(k) == 1) then
                                        interior_points = interior_points + 1
                                else if(point%flag_1(k) == 2) then
                                        outer_points = outer_points + 1
                                end if

                                if(point%flag_2(k) > shapes)then
                                        print*,"shapes value wrong, check again"
                                        stop
                                end if
                                         
        
                                
                        enddo
                else
                        do k = 1, max_points
        
                                read(101,*) point%local_id(k),point%global_id(k),point%x(k),&
                                & point%y(k),point%nx(k),point%ny(k), point%flag_1(k),point%flag_2(k),point%nbhs(k),&
                                & (point%conn(k,r),r=1,point%nbhs(k))
                                
                        !Storing the count for the point types
                                if(point%flag_1(k) == 1) then
                                        wall_points = wall_points + 1
                                else if(point%flag_1(k) == 2) then
                                        interior_points = interior_points + 1
                                else if(point%flag_1(k) == 3) then
                                        outer_points = outer_points + 1
                                end if
        
                                
                        enddo
                end if


                allocate(wall_points_index(wall_points))
                allocate(interior_points_index(interior_points))
                allocate(outer_points_index(outer_points))

                if(old_format == 0)then
                        wall_temp = 0
                        interior_temp = 0
                        outer_temp = 0
                        !Storing indices of the point definitions
                        do k = 1,max_points
                                if(point%flag_1(k) == 0) then
                                        wall_temp = wall_temp+1
                                        wall_points_index(wall_temp) = k
                                else if(point%flag_1(k) == 1) then
                                        interior_temp = interior_temp+1 
                                        interior_points_index(interior_temp) = k
                                else if(point%flag_1(k) == 2) then
                                        outer_temp = outer_temp+1
                                        outer_points_index(outer_temp) = k
                                end if
                        end do
                else
                        wall_temp = 0
                        interior_temp = 0
                        outer_temp = 0
                        !Storing indices of the point definitions
                        do k = 1,max_points
                                if(point%flag_1(k) == 1) then
                                        wall_temp = wall_temp+1
                                        wall_points_index(wall_temp) = k
                                else if(point%flag_1(k) == 2) then
                                        interior_temp = interior_temp+1 
                                        interior_points_index(interior_temp) = k
                                else if(point%flag_1(k) == 3) then
                                        outer_temp = outer_temp+1
                                        outer_points_index(outer_temp) = k
                                end if
                        end do
                end if

                CLOSE(UNIT=101)

        end subroutine 

        subroutine dealloc_points()
                implicit none

                
                deallocate(point%x)
                deallocate(point%y)
                deallocate(point%local_id)
                deallocate(point%global_id)
                deallocate(point%flag_1)
                deallocate(point%flag_2)
                deallocate(point%nbhs)
                deallocate(point%conn)
                deallocate(point%nx)
                deallocate(point%ny)
                deallocate(point%left)
                deallocate(point%right)

                deallocate(wall_points_index)
                deallocate(interior_points_index)
                deallocate(outer_points_index)

        end subroutine

end module 
