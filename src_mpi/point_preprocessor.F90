module point_preprocessor_mod

        use data_structure_mod

contains

        subroutine read_input_point_data()

#include <petsc/finclude/petscsys.h>

                use petscsys

                implicit none

                integer:: i, k, r
                integer :: wall_temp,outer_temp,interior_temp
                character(len=64) :: part_grid
                !TODO : write a script to convert integer to string up to 4 digits
                character(len=10) :: itos

                if (rank==0) print*,'Reading points'

                part_grid = 'partGrid0'
                if (proc>1) part_grid = 'par/partGrid'//trim(itos(1,rank))

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
                             !   wall_points_index(wall_points) = k
                        else if(point(k)%flag_1 == 2) then
                                interior_points = interior_points + 1
                              !  interior_points_index(interior_points) = k
                        else if(point(k)%flag_1 == 3) then
                                outer_points = outer_points + 1
                               ! outer_points_index(outer_points) = k
                        end if

                        

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

                        

                if (proc > 1) then
                        allocate(pghost(ghost_points))
                        do k=local_points+1,max_points
                                read(101,*) point(k)%local_id,pghost(k-local_points),&
                                & point(k)%x,point(k)%y
                        enddo
                end if


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
