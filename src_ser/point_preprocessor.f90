module point_preprocessor_mod

        use data_structure_mod

contains

        subroutine read_input_point_data()



                implicit none

                integer:: i, k, r
                integer :: wall_temp,outer_temp,interior_temp
                character(len=64) :: part_grid
                !TODO : write a script to convert integer to string up to 4 digits
                character(len=10) :: itos


                part_grid = 'partGrid0'

		OPEN(UNIT=101,FILE=trim(part_grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

                !TODO: Add asserts

!               input file format : 
!               max_points local_points ghost_points
!               local_id global_id x_cord y_cord flag_1 flag_2 num_Connectivity connectivity[]

		read(101,*) max_points, local_points, ghost_points
                allocate(point(max_points))

		wall_points = 0
		interior_points = 0
		outer_points = 0
                shape_points = 0
                
		do k = 1, local_points

                        read(101,*) point(k)%local_id,point(k)%global_id,point(k)%x,&
                        & point(k)%y, point(k)%flag_1,point(k)%flag_2,point(k)%nbhs,&
                        & (point(k)%conn(r),r=1,point(k)%nbhs)

                        if(point(k)%flag_1 == 1) then
                                wall_points = wall_points + 1
                        else if(point(k)%flag_1 == 2) then
                                interior_points = interior_points + 1
                        else if(point(k)%flag_1 == 3) then
                                outer_points = outer_points + 1
                        end if

                        
!TODO: They will be greater than zero. 
                        if(point(k)%flag_2 > 0) then
                                shape_points(point(k)%flag_2) = shape_points(point(k)%flag_2)+1
                                shape_points_index(point(k)%flag_2,shape_points(point(k)%flag_2))=k
                        end if

		enddo
                allocate(wall_points_index(wall_points))
                allocate(interior_points_index(interior_points))
                allocate(outer_points_index(outer_points))

                wall_temp = 0
                interior_temp = 0
                outer_temp = 0
                do k = 1,local_points
                        if(point(k)%flag_1 == 1) then
                                wall_temp = wall_temp+1
                                wall_points_index(wall_temp) = k
                        else if(point(k)%flag_1 == 2) then
                                interior_temp = interior_temp+1 
                                interior_points_index(interior_temp) = k
                        else if(point(k)%flag_1 == 3) then
                                outer_temp = outer_temp+1
                                outer_points_index(outer_temp) = k
                        end if
                end do

                        



!	The above lines of the code will remain the same for all test cases. 
!	However, depending on the number of shapes in the geometry, we divide
!	the wall points into respective shape points ..	
!		do r = 1, shapes
!			do i = 1, shape_points(r)
!				read(101,*) shape_points_index(r, i) 				
!			enddo
!		enddo

		CLOSE(UNIT=101)

	end subroutine 
				
end module 		
