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
                character(len=10) :: itos

                if (rank==0) print*,'Reading points'
!TODO Make preproc to name ser grid as partgrid and par as partgrid0....
                part_grid = 'partGrid0'
                if (proc>1) part_grid = 'par/partGrid'//trim(itos(1,rank))

                OPEN(UNIT=101,FILE=trim(part_grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

                !TODO: Add asserts

!               input file format : 
!               max_points local_points ghost_points
!               local_id global_id x_cord y_cord flag_1 flag_2 num_Connectivity connectivity[]

                read(101,*) max_points, local_points, ghost_points
                allocate(p%x(max_points))
                allocate(p%y(max_points))
                allocate(p%local_id(max_points))
                allocate(p%global_id(max_points))
                allocate(p%flag_1(max_points))
                allocate(p%flag_2(max_points))
                allocate(p%nbhs(max_points))
                allocate(p%conn(max_points,15))

                wall_points = 0
                interior_points = 0
                outer_points = 0
                shape_points = 0
                
                do k = 1, local_points

                        read(101,*) p%local_id(k),p%global_id(k),p%x(k),&
                        & p%y(k), p%flag_1(k),p%flag_2(k),p%nbhs(k),&
                        & (p%conn(k,r),r=1,p%nbhs(k))
                        
                !Storing the count for the point types
                        if(p%flag_1(k) == 1) then
                                wall_points = wall_points + 1
                        else if(p%flag_1(k) == 2) then
                                interior_points = interior_points + 1
                        else if(p%flag_1(k) == 3) then
                                outer_points = outer_points + 1
                        end if

                        
                !Storing shape point indices for all the shape types
                        if(p%flag_2(k) > 0) then
                                shape_points(p%flag_2(k)) = shape_points(p%flag_2(k))+1
                                shape_points_index(p%flag_2(k),shape_points(p%flag_2(k)))=k
                        end if

                enddo


                allocate(wall_points_index(wall_points))
                allocate(interior_points_index(interior_points))
                allocate(outer_points_index(outer_points))

                wall_temp = 0
                interior_temp = 0
                outer_temp = 0
                !Storing indices of the point definitions
                do k = 1,local_points
                        if(p%flag_1(k) == 1) then
                                wall_temp = wall_temp+1
                                wall_points_index(wall_temp) = k
                        else if(p%flag_1(k) == 2) then
                                interior_temp = interior_temp+1 
                                interior_points_index(interior_temp) = k
                        else if(p%flag_1(k) == 3) then
                                outer_temp = outer_temp+1
                                outer_points_index(outer_temp) = k
                        end if
                end do

                        

                if (proc > 1) then
                        allocate(pghost(ghost_points))

                        do k=local_points+1,max_points
                                read(101,*) p%local_id(k),pghost(k-local_points),&
                                & p%x(k),p%y(k),p%flag_1(k),p%flag_2(k)

                        end do
                end if









                CLOSE(UNIT=101)

        end subroutine 

end module 
