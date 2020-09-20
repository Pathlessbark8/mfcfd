module point_preprocessor_mod

    use data_structure_mod
    use petsc_data_structure_mod
    USE ReadH5dataset

contains

    subroutine read_input_point_data()

#include <petsc/finclude/petscsys.h>

        use petscsys

        implicit none

        integer:: i, k, r, nproc
        integer :: wall_temp,outer_temp,interior_temp,shape_temp
        character(len=64) :: part_grid
        character(len=10) :: itos

        part_grid = 'point/partGrid'
        if (proc>1) part_grid = 'point/partGrid'//trim(itos(4,rank))

        OPEN(UNIT=101,FILE=trim(part_grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

        !TODO: Add asserts

        read(101,*) nproc, max_points, local_points, ghost_points
        if(proc .ne. nproc) then
            SETERRA(PETSC_COMM_WORLD,1,'check number of partitions and proc')
        end if
        allocate(point%x(max_points))
        allocate(point%y(max_points))
        allocate(point%flag_1(max_points))
        allocate(point%flag_2(max_points))
        allocate(point%nx(max_points))
        allocate(point%ny(max_points))
        allocate(point%qtdepth(max_points))
        allocate(point%nbhs(max_points))
        allocate(point%conn(max_points,20))
        allocate(point%min_dist(max_points))
        allocate(point%left(max_points))
        allocate(point%right(max_points))
        allocate(point%original_id(local_points))

        wall_points = 0
        interior_points = 0
        outer_points = 0
        shape_points = 0

        if(format == 1) then


            do k = 1, local_points
            
                read(101,*)point%original_id(k), point%x(k), point%y(k), &
                & point%left(k),point%right(k), point%flag_1(k),point%flag_2(k), &
                & point%min_dist(k), point%nbhs(k), (point%conn(k,r),r=1,point%nbhs(k))

            !Storing the count for the point types
                if(point%flag_1(k) == 0) then
                    wall_points = wall_points + 1
                else if(point%flag_1(k) == 1) then
                    interior_points = interior_points + 1
                else if(point%flag_1(k) == 2) then
                    outer_points = outer_points + 1
                end if

                if(point%flag_2(k) > 0) then
                    if(point%flag_2(k) > shapes)then
                        SETERRA(PETSC_COMM_WORLD,1,'shapes value wrong, check again')
                    end if
                    shape_points = shape_points + 1
                end if

            enddo
        else if(format == 2) then ! quadtree format
            do k = 1, local_points
            
                read(101,*)point%original_id(k), point%x(k), point%y(k), &
                & point%left(k),point%right(k), point%flag_1(k),point%flag_2(k), &
                & point%nx(k), point%ny(k), point%qtdepth(k),point%min_dist(k), &
                & point%nbhs(k), (point%conn(k,r),r=1,point%nbhs(k))
                
            !Storing the count for the point types
                if(point%flag_1(k) == 0) then
                    wall_points = wall_points + 1
                else if(point%flag_1(k) == 1) then
                    interior_points = interior_points + 1
                else if(point%flag_1(k) == 2) then
                    outer_points = outer_points + 1
                end if
            
                if(point%flag_2(k) > 0) then
                      if(point%flag_2(k) > shapes)then
                          print*,"shapes value wrong, check again"
                          stop
                      end if
                      shape_points = shape_points + 1
                end if
            
            enddo

        end if

        allocate(wall_points_index(wall_points))
        allocate(interior_points_index(interior_points))
        allocate(outer_points_index(outer_points))
        allocate(shape_points_index(shape_points))


        wall_temp = 0
        interior_temp = 0
        outer_temp = 0
        shape_temp = 0
        !Storing indices of the point definitions
        do k = 1,local_points
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

            if(point%flag_2(k) > 0) then
                shape_temp = shape_temp+1
                shape_points_index(shape_temp) = k
            end if
        end do

        if (proc > 1) then
            allocate(pghost(ghost_points))

            do k= 1, ghost_points
                read(101,*) pghost(k),&
                & point%x(local_points + k),point%y(local_points + k), &
                & point%min_dist(local_points + k)

            end do
        end if

        CLOSE(UNIT=101)

    end subroutine 

    subroutine dealloc_points()
        implicit none

        deallocate(point%x)
        deallocate(point%y)
        deallocate(point%flag_1)
        deallocate(point%flag_2)
        deallocate(point%nx)
        deallocate(point%ny)
        deallocate(point%nbhs)
        deallocate(point%qtdepth)
        deallocate(point%conn)
        deallocate(point%min_dist)
        deallocate(point%left)
        deallocate(point%right)
        deallocate(point%original_id)

        deallocate(wall_points_index)
        deallocate(interior_points_index)
        deallocate(outer_points_index)
        
        if(allocated(pghost)) deallocate(pghost)

    end subroutine

    subroutine read_hdf5input_point_data()

#include <petsc/finclude/petscsys.h>

        use petscsys

        implicit none

        integer:: i, k, r, nproc
        integer :: wall_temp,outer_temp,interior_temp,shape_temp
        character(len=65) :: part_grid, dataset_string

        integer, dimension(:),pointer :: nbh_array
        real*8 :: nx

        character(len=10) :: itos, itos_unpad
        character(len=10) :: main_group, total_attribute, ghost_attribute, local_attribute, point_string, rank_string 
        character(len=10) :: flag1_s, flag2_s, left_s, right_s, min_dist_s, nx_s
        character(len=10) :: ny_s, qtdepth_s, x_s, y_s
        character(len=15) :: nbhs_counts_s

        part_grid = 'point/point.h5'
        if (proc>1) part_grid = 'point/partGrid'//trim(itos(4,rank))

        main_group = '/'//trim(itos_unpad(rank+1))
        ghost_attribute = '/ghost'
        local_attribute = '/local'
        total_attribute = '/total'
        flag1_s = '/flag1'  
        flag2_s = '/flag2'  
        left_s = '/left'  
        right_s = '/right'  
        min_dist_s = '/min_dist'  
        nbhs_counts_s = '/nbhs_count'  
        nx_s = '/nx' 
        ny_s = '/ny'  
        qtdepth_s = '/qtdepth'  
        x_s = '/x'  
        y_s = '/y' 

        CALL H5Read_initfile(part_grid, main_group)

        ! CALL H5ReadAttribute(part_grid, main_group//total_attribute, max_points)
        ! CALL H5ReadAttribute(part_grid, main_group//total_attribute, max_points)
        ! CALL H5ReadAttribute(part_grid, "/1/ghost", ghost_points)
        ! CALL H5ReadAttribute(part_grid, main_group//local_attribute, local_points)

        local_points = max_points
        write(*,*) local_points, max_points

        allocate(point%x(max_points))
        allocate(point%y(max_points))
        allocate(point%flag_1(max_points))
        allocate(point%flag_2(max_points))
        allocate(point%nx(max_points))
        allocate(point%ny(max_points))
        allocate(point%qtdepth(max_points))
        allocate(point%nbhs(max_points))
        allocate(point%conn(max_points,20))
        allocate(point%min_dist(max_points))
        allocate(point%left(max_points))
        allocate(point%right(max_points))
        allocate(point%original_id(local_points))

        wall_points = 0
        interior_points = 0
        outer_points = 0
        shape_points = 0

        ! do k = 1, 100
            k = 1
            point_string = '/'//trim(itos_unpad(k))
            dataset_string = main_group//local_attribute//point_string

            point%original_id(k) = k
            ! CALL H5ReadAttribute(part_grid, dataset_string//x_s, point%x(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//y_s, point%y(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//left_s, point%left(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//right_s, point%right(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//flag1_s, point%flag_1(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//flag2_s, point%flag_2(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//min_dist_s, point%min_dist(k))
            CALL H5ReadAttribute(part_grid, dataset_string//nbhs_counts_s, point%nbhs(k))
            ! CALL H5ReadAttribute(part_grid, dataset_string//qtdepth_s, point%qtdepth(k))
            CALL H5ReadDataset(part_grid, dataset_string, nbh_array)

            do r=1,point%nbhs(k)
                WRITE(*,*) nbh_array(r)
            end do
            ! !Storing the count for the point types
            ! if(point%flag_1(k) == 0) then
            !     wall_points = wall_points + 1
            ! else if(point%flag_1(k) == 1) then
            !     interior_points = interior_points + 1
            ! else if(point%flag_1(k) == 2) then
            !     outer_points = outer_points + 1
            ! end if

            ! if(point%flag_2(k) > 0) then
            !     if(point%flag_2(k) > shapes)then
            !         SETERRA(PETSC_COMM_WORLD,1,'shapes value wrong, check again')
            !     end if
            !     shape_points = shape_points + 1
            ! end if

        ! end do    

        allocate(wall_points_index(wall_points))
        allocate(interior_points_index(interior_points))
        allocate(outer_points_index(outer_points))
        allocate(shape_points_index(shape_points))


        ! wall_temp = 0
        ! interior_temp = 0
        ! outer_temp = 0
        ! shape_temp = 0
        ! !Storing indices of the point definitions
        ! do k = 1,local_points
        !     if(point%flag_1(k) == 0) then
        !         wall_temp = wall_temp+1
        !         wall_points_index(wall_temp) = k
        !     else if(point%flag_1(k) == 1) then
        !         interior_temp = interior_temp+1 
        !         interior_points_index(interior_temp) = k
        !     else if(point%flag_1(k) == 2) then
        !         outer_temp = outer_temp+1
        !         outer_points_index(outer_temp) = k
        !     end if

        !     if(point%flag_2(k) > 0) then
        !         shape_temp = shape_temp+1
        !         shape_points_index(shape_temp) = k
        !     end if
        ! end do  

        CALL H5Read_closefile()
    end subroutine 

end module 
