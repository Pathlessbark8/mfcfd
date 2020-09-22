module point_preprocessor_mod

        use data_structure_mod
        use device_data_structure_mod
        USE HDF5

contains

        subroutine read_input_point_data()

                implicit none

                integer:: i, k, r
                integer :: wall_temp,outer_temp,interior_temp,shape_temp
                character(len=64) :: grid
                real*8 :: dummy

                grid = 'point/partGrid'

                OPEN(UNIT=101,FILE=trim(grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

                read(101,*) max_points

                allocate(point%x(max_points))
                allocate(point%y(max_points))
                allocate(point%flag_1(max_points))
                allocate(point%flag_2(max_points))
                allocate(point%min_dist(max_points))
                allocate(point%nbhs(max_points))
                allocate(point%conn(max_points,25))
                allocate(point%nx(max_points))
                allocate(point%ny(max_points))
                allocate(point%left(max_points))
                allocate(point%right(max_points))
                allocate(point%qtdepth(max_points))
                
                ! device variables
                allocate(point_d%x(2,max_points))
                allocate(point_d%flag(max_points))
                allocate(point_d%nbhs(max_points))
                allocate(point_d%conn(max_points,25))
                allocate(point_d%nx(2,max_points))
                allocate(point_d%min_dist(max_points))
                
                wall_points = 0
                shape_points = 0
                interior_points = 0
                outer_points = 0
                ! legacy format
                if (file_format == 1) then
                        do k = 1, max_points
        
                                read(101,*) point%x(k),&
                                & point%y(k),point%left(k),point%right(k), &
                                & point%flag_1(k),point%flag_2(k),point%min_dist(k), &
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
                ! new format from quadtree code
                elseif (file_format == 2) then
                        do k = 1, max_points
        
                                read(101,*) point%x(k),&
                                & point%y(k),point%left(k),point%right(k),point%flag_1(k), &
                                & point%flag_2(k), point%nx(k), point%ny(k), &
                                & point%qtdepth(k), point%min_dist(k), point%nbhs(k),&
                                & (point%conn(k,r),r=1,point%nbhs(k))
                                
                        !Storing the count for the point types
                                if(point%flag_1(k) == 0) then
                                        wall_points = wall_points + 1
                                else if(point%flag_1(k) == 1) then
                                        interior_points = interior_points + 1
                                else if(point%flag_1(k) == 2) then
                                        outer_points = outer_points + 1
                                end if

                                if(point%nbhs(k) == 0) then
                                        print*,"No neighbours for point:", k
                                        stop
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
                allocate(shape_points_index(shape_points))
                allocate(interior_points_index(interior_points))
                allocate(outer_points_index(outer_points))

                wall_temp = 0
                shape_temp = 0
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
                        
                        if(point%flag_2(k) > 0) then
                                shape_temp = shape_temp+1
                                shape_points_index(shape_temp) = k
                        end if
                end do

                CLOSE(UNIT=101)

        end subroutine 

        subroutine dealloc_points()
                implicit none

                
                deallocate(point%x)
                deallocate(point%y)
                deallocate(point%flag_1)
                deallocate(point%flag_2)
                deallocate(point%min_dist)
                deallocate(point%nbhs)
                deallocate(point%conn)
                deallocate(point%nx)
                deallocate(point%ny)
                deallocate(point%left)
                deallocate(point%right)
                deallocate(point%qtdepth)

                deallocate(wall_points_index)
                deallocate(shape_points_index)
                deallocate(interior_points_index)
                deallocate(outer_points_index)

        end subroutine

    subroutine read_hdf5input_point_data()

        implicit none

        integer:: i, k, r
        integer :: wall_temp,outer_temp,interior_temp,shape_temp
        character(len=65) :: part_grid, dataset_string

        CHARACTER(LEN = 1)  :: rootname
        character(len=10) :: itos, itos_unpad
        character(len=10) :: main_group, total_attribute, ghost_attribute, local_attribute, point_string, rank_string 
        character(len=10) :: val1_s, val2_s
        character(len=15) :: local_group, ghost_group
        CHARACTER(LEN=100)  :: ErrorMessage

        INTEGER(HID_T)      :: file_id
        INTEGER(HID_T)      :: root_id
        INTEGER(HID_T)      :: group_id, localgroup_id, ghostgroup_id
        INTEGER(HID_T)      :: p_id
        INTEGER(HID_T)      :: d_id, a_id
        INTEGER(HID_T)      :: type_id
        INTEGER(HID_T)      :: dspace
        INTEGER             :: ErrorFlag
        INTEGER             :: AllocStat
        INTEGER             :: H5dataset
        INTEGER(hsize_t), DIMENSION(1)                       :: dims
        INTEGER, DIMENSION(:), ALLOCATABLE, TARGET          :: H51DIntegerdataset
        INTEGER, dimension(:), pointer :: nbh_array 
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, TARGET     :: H51DDoubledataset

        part_grid = 'point/point.h5'
        ! if (proc>1) part_grid = 'point/partGrid'//trim(itos(4,rank))

        main_group = '/'//trim(itos_unpad(rank+1))
        ghost_attribute = 'ghost'
        local_attribute = 'local'
        total_attribute = 'total'
        val1_s = 'val1'
        val2_s = 'val2'

        CALL h5open_f(ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error initialising HDF routines"
            return
        ENDIF
        
        ! CALL DebugMessage(" --- Opening HDF file")
        
        CALL h5fopen_f(part_grid, H5F_ACC_RDONLY_F, file_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening HDF file"
            return
        ENDIF
        
        ! Open the root group.
        rootname = "/"
        CALL h5gopen_f(file_id,rootname,root_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening root group"
            return
        ENDIF

        ! Main group
        CALL h5gopen_f(root_id, main_group, group_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening main group "
            return
        ENDIF

        CALL h5aopen_name_f(group_id, total_attribute, a_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening total attribute "
            return
        ENDIF
        dims=1
        CALL h5aread_f(a_id, H5T_NATIVE_INTEGER, max_points, dims, ErrorFlag)
        CALL h5aclose_f(a_id,ErrorFlag)


        CALL h5aopen_name_f(group_id, local_attribute, a_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening local attribute "
            return
        ENDIF
        dims=1
        CALL h5aread_f(a_id, H5T_NATIVE_INTEGER, local_points, dims, ErrorFlag)
        CALL h5aclose_f(a_id,ErrorFlag)

        CALL h5aopen_name_f(group_id, ghost_attribute, a_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening ghost attribute "
            return
        ENDIF
        dims=1
        CALL h5aread_f(a_id, H5T_NATIVE_INTEGER, ghost_points, dims, ErrorFlag)
        CALL h5aclose_f(a_id,ErrorFlag)

        allocate(point%x(max_points))
        allocate(point%y(max_points))
        allocate(point%flag_1(max_points))
        allocate(point%flag_2(max_points))
        allocate(point%nx(max_points))
        allocate(point%ny(max_points))
        allocate(point%qtdepth(max_points))
        allocate(point%nbhs(max_points))
        allocate(point%conn(max_points,25))
        allocate(point%min_dist(max_points))
        allocate(point%left(max_points))
        allocate(point%right(max_points))
                        ! device variables
        allocate(point_d%x(2,max_points))
        allocate(point_d%flag(max_points))
        allocate(point_d%nbhs(max_points))
        allocate(point_d%conn(max_points,25))
        allocate(point_d%nx(2,max_points))
        allocate(point_d%min_dist(max_points))

        ALLOCATE(H51DIntegerdataset(7))
        ALLOCATE(H51DDoubledataset(5))
        ALLOCATE(nbh_array(20))

        wall_points = 0
        interior_points = 0
        outer_points = 0
        shape_points = 0

        local_group = trim(main_group)//'/local'
        CALL h5gopen_f(group_id, local_group, localgroup_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening local group "
            return
        ENDIF

        do k = 1, max_points
            ! k = 1
            point_string = '/'//trim(itos_unpad(k))
            dataset_string = trim(local_group)//point_string

            CALL h5dopen_f(localgroup_id , dataset_string , d_id, ErrorFlag)
            IF (ErrorFlag.lt.0) THEN
                ErrorMessage=" *** Error opening dataset "
                return
            ENDIF
            dims=20
            CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, nbh_array, dims, ErrorFlag)
            ! WRITE(*,*) nbh_array

            CALL h5aopen_name_f(d_id, val2_s, a_id, ErrorFlag)
            IF (ErrorFlag.lt.0) THEN
                ErrorMessage=" *** Error opening attribute "
                return
            ENDIF
            dims=7
            CALL h5aread_f(a_id, H5T_NATIVE_INTEGER, H51DIntegerdataset, dims, ErrorFlag)
            CALL h5aclose_f(a_id,ErrorFlag)

            ! WRITE(*,*) H51DIntegerdataset

            CALL h5aopen_name_f(d_id, val1_s, a_id, ErrorFlag)
            IF (ErrorFlag.lt.0) THEN
                ErrorMessage=" *** Error opening attribute "
                return
            ENDIF
            dims=5
            CALL h5aread_f(a_id, H5T_NATIVE_DOUBLE, H51DDoubledataset, dims, ErrorFlag)
            CALL h5aclose_f(a_id,ErrorFlag)

            ! WRITE(*,*) H51DDoubledataset

            CALL h5dclose_f(d_id,ErrorFlag)

            point%x(k) = H51DDoubledataset(1)
            point%y(k) = H51DDoubledataset(2)
            point%left(k) = H51DIntegerdataset(2)
            point%right(k) = H51DIntegerdataset(3)
            point%flag_1(k) = H51DIntegerdataset(5)
            point%flag_2(k) = H51DIntegerdataset(6)
            point%nx(k) = H51DDoubledataset(3)
            point%ny(k) = H51DDoubledataset(4)
            point%qtdepth(k) = H51DIntegerdataset(4)
            point%min_dist(k) = H51DDoubledataset(5)
            point%nbhs(k) = H51DIntegerdataset(7)
            do r=1, point%nbhs(k)
                point%conn(k,r) = nbh_array(r)
            end do

            ! !Storing the count for the point types
            if(point%flag_1(k) == 0) then
                wall_points = wall_points + 1
            else if(point%flag_1(k) == 1) then
                interior_points = interior_points + 1
            else if(point%flag_1(k) == 2) then
                outer_points = outer_points + 1
            end if

            if(point%nbhs(k) == 0) then
                print*,"No neighbours for point:", k
                stop
            end if

            if(point%flag_2(k) > 0) then
                if(point%flag_2(k) > shapes)then
                    print*,"shapes value wrong, check again"
                    stop
                end if
                shape_points = shape_points + 1
            end if

        end do    

        CALL h5gclose_f(localgroup_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error closing local group"
            return
        ENDIF

        allocate(wall_points_index(wall_points))
        allocate(interior_points_index(interior_points))
        allocate(outer_points_index(outer_points))
        allocate(shape_points_index(shape_points))

        deallocate(nbh_array)
        deallocate(H51DIntegerdataset)
        deallocate(H51DDoubledataset)


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

        CALL h5gclose_f(group_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error closing main group"
            return
        ENDIF
        

        CALL h5gclose_f(root_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error closing root group"
            return
        ENDIF
        
        ! Close the HDF file and HDF interface
        
        ! CALL DebugMessage(" --- Closing HDF file")
        CALL h5fclose_f(file_id,ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error closing HDF file"
            return
        ENDIF
        
        CALL h5close_f(ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error closing HDF routines"
            return
        ENDIF
    end subroutine 


end module 
