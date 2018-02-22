module point_preprocessor_mod
!
	use data_structure_mod
!
	contains

		subroutine read_input_point_data()
!
!
		implicit none

		integer i, k, r

!		NOTE : change file path accordingly to single partition file		
		OPEN(UNIT=101,FILE="./../../metisPre/mfapre/partGrid0",FORM="FORMATTED",STATUS="OLD",ACTION="READ")

		! input file format : 
        ! max_points local_points
        ! local_id global_id x_cord y_cord flag_1 flag_2 num_Connectivity connectivity[]

		read(101,*) max_points, local_points

		wall_points = 0
		interior_points = 0
		outer_points = 0
		do i = 1, max_points
				read(101,*) k, point(k)%global_id, point(k)%x, point(k)%y, &
				& point(k)%flag_1, point(k)%flag_2, point(k)%nbhs, (point(k)%conn(r),r=1,point(k)%nbhs)
				IF(point(k)%flag_1 == 1) THEN
					wall_points = wall_points + 1
				ELSE IF(point(k)%flag_1 == 2) THEN
					interior_points = interior_points + 1
				ELSE IF(point(k)%flag_1 == 3) THEN
					outer_points = outer_points + 1
				END IF
		enddo		

!	The above lines of the code will remain the same for all test cases. 
!	However, depending on the number of shapes in the geometry, we divide
!	the wall points into respective shape points ..
!	
!	
		! do r = 1, shapes
		! 	do i = 1, shape_points(r)
		! 		read(101,*) shape_points_index(r, i) 				
		! 	enddo
		! enddo
!
!		
!
		CLOSE(UNIT=101)
!
!
	end subroutine 
!
!					
end module 		
