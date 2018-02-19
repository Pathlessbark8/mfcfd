module point_preprocessor_mod
!
	use data_structure_mod
!
	contains
	!
	!
!
!
!
	subroutine read_input_point_data()
!
!
		implicit none
!
		integer i, k, r
!
!
!		OPEN(UNIT=101,FILE="./unstructured-grid-4733/input-file",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
!		OPEN(UNIT=101,FILE="./structured-grid-640-240/input-file",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
!		OPEN(UNIT=101,FILE="./structured-grid-320-120/input-file",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
		OPEN(UNIT=101,FILE="./structured-grid-160-60/input-file",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
!		OPEN(UNIT=102,FILE="check-input",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
!		
!
!		read(101,*) max_points, interior_points, outer_points, wall_points, shapes, (shape_points(r), r=1, shapes)
		!
		do i = 1, max_points
				read(101,*) k, point(k)%x, point(k)%y, point(k)%nbhs, (point(k)%conn(r),r=1,point(k)%nbhs)
		enddo		
!
		do i = 1, interior_points
			read(101,*) interior_points_index(i) 
		enddo		
!
		do i = 1, outer_points
			read(101,*) outer_points_index(i)
		enddo	
!
		do i = 1, wall_points
			read(101,*) wall_points_index(i) 
		enddo
		!
!
!	The above lines of the code will remain the same for all test cases. 
!	However, depending on the number of shapes in the geometry, we divide
!	the wall points into respective shape points ..
!	
!	
		do r = 1, shapes
			do i = 1, shape_points(r)
				read(101,*) shape_points_index(r, i) 				
			enddo
		enddo
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
