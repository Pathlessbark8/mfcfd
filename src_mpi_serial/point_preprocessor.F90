module point_preprocessor_mod

    use data_structure_mod
    use petsc_data_structure_mod
    USE HDF5
    ! USE ReadH5dataset

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

        CHARACTER(LEN = 1)  :: rootname
        character(len=10) :: itos, itos_unpad
        character(len=10) :: main_group, total_attribute, ghost_attribute, local_attribute, point_string, rank_string 
        character(len=10) :: val1_s, val2_s
        character(len=15) :: local_group, ghost_group, local_dataset, ghost_dataset
        CHARACTER(LEN=100)  :: ErrorMessage

        INTEGER(HID_T)      :: file_id
        INTEGER(HID_T)      :: root_id
        INTEGER(HID_T)      :: group_id, localgroup_id, ghostgroup_id, localdatset_id, ghostdataset_id
        INTEGER(HID_T)      :: p_id
        INTEGER(HID_T)      :: d_id, a_id
        INTEGER(HID_T)      :: type_id
        INTEGER(HID_T)      :: dspace
        INTEGER             :: ErrorFlag, ndims
        INTEGER             :: AllocStat
        INTEGER             :: H5dataset
        INTEGER             :: maxdims = 4
        INTEGER(hsize_t), DIMENSION(1)                 :: dims
        INTEGER(hsize_t), DIMENSION(:), ALLOCATABLE    :: main_dims
        INTEGER(hsize_t), DIMENSION(4)                 :: datadims
        INTEGER(hsize_t), DIMENSION(4)                 :: maxdatadims
        INTEGER, DIMENSION(:), ALLOCATABLE, TARGET          :: H51DIntegerdataset
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, TARGET     :: H51DDoubledataset
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, TARGET     :: H52DDoubledataset

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
        allocate(point%conn(max_points,20))
        allocate(point%min_dist(max_points))
        allocate(point%left(max_points))
        allocate(point%right(max_points))
        allocate(point%original_id(local_points))

        wall_points = 0
        interior_points = 0
        outer_points = 0
        shape_points = 0

        local_dataset = trim(main_group)//'/local'

        CALL h5dopen_f(group_id , local_dataset , d_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening dataset "
            return
        ENDIF

        ! Find space of 2D array

        CALL h5dget_space_f(d_id, dspace, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error determining dataspace"
            return
        ENDIF
     
        CALL h5sget_simple_extent_ndims_f(dspace, ndims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
           ErrorMessage=" *** Error determining dataspace dimensionality"
           return
        ENDIF

        ALLOCATE(main_dims(ndims),stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error allocating dims"
           return
        ENDIF

        CALL h5sget_simple_extent_dims_f(dspace, datadims, maxdatadims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
           ErrorMessage=" *** Error determining dataspace size"
           return
        ENDIF
    
        main_dims = datadims
    
        ! Read the dataset and transfer the result
    
        ALLOCATE(H52DDoubledataset(main_dims(1), main_dims(2)), stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error allocating H5dataset"
           return
        ENDIF

        CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H52DDoubledataset, dims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error reading data"
            return
        ENDIF

        do k = 1, local_points
        !     ! k = 1
        !     ! WRITE(*,*) nbh_array

            point%original_id(k) = int(H52DDoubledataset(1,k))
            point%x(k) = H52DDoubledataset(2,k)
            point%y(k) = H52DDoubledataset(3,k)
            point%left(k) = H52DDoubledataset(7,k)
            point%right(k) = H52DDoubledataset(8,k)
            point%flag_1(k) = H52DDoubledataset(10,k)
            point%flag_2(k) = H52DDoubledataset(11,k)
            point%nx(k) = H52DDoubledataset(4,k)
            point%ny(k) = H52DDoubledataset(5,k)
            point%qtdepth(k) = H52DDoubledataset(9,k)
            point%min_dist(k) = H52DDoubledataset(6,k)
            point%nbhs(k) = H52DDoubledataset(12,k)
            do r=1, point%nbhs(k)
                point%conn(k,r) = H52DDoubledataset(12+r,k)
            end do

            ! !Storing the count for the point types
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

        end do    

        DEALLOCATE(H52DDoubledataset,stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error deallocating H52DDoubledatasets"
           return
        ENDIF    

        DEALLOCATE(main_dims, stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error deallocating dims"
           return
        ENDIF

        CALL h5dclose_f(d_id, ErrorFlag)

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

            if(point%original_id(k)== 9521) then
                WRITE(*,*) point%x(k-1), point%x(k+1)
            endif
        end do  

        ghost_dataset = trim(main_group)//'/ghost'

        CALL h5dopen_f(group_id , ghost_dataset , d_id, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error opening dataset "
            return
        ENDIF

        ! Find space of 2D array

        CALL h5dget_space_f(d_id, dspace, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error determining dataspace"
            return
        ENDIF
     
        CALL h5sget_simple_extent_ndims_f(dspace, ndims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
           ErrorMessage=" *** Error determining dataspace dimensionality"
           return
        ENDIF

        ALLOCATE(main_dims(ndims),stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error allocating dims"
           return
        ENDIF

        CALL h5sget_simple_extent_dims_f(dspace, datadims, maxdatadims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
           ErrorMessage=" *** Error determining dataspace size"
           return
        ENDIF
    
        main_dims = datadims
    
        ! Read the dataset and transfer the result
    
        ALLOCATE(H52DDoubledataset(main_dims(1), main_dims(2)), stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error allocating H5dataset"
           return
        ENDIF

        CALL h5dread_f(d_id, H5T_NATIVE_DOUBLE, H52DDoubledataset, dims, ErrorFlag)
        IF (ErrorFlag.lt.0) THEN
            ErrorMessage=" *** Error reading data"
            return
        ENDIF

        if (proc > 1) then
            allocate(pghost(ghost_points))

            do k = 1, ghost_points
    
                pghost(k) = int(H52DDoubledataset(1,k))

                point%x(local_points + k) = H52DDoubledataset(2,k)
                point%y(local_points + k) = H52DDoubledataset(3,k)
                point%min_dist(local_points + k) = H52DDoubledataset(4,k)
            end do

        end if

        DEALLOCATE(H52DDoubledataset,stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error deallocating H52DDoubledatasets"
           return
        ENDIF    

        DEALLOCATE(main_dims, stat=AllocStat)
        IF ( AllocStat.ne.0 ) THEN
           ErrorFlag=-1
           ErrorMessage=" *** Error deallocating dims"
           return
        ENDIF

        CALL h5dclose_f(d_id, ErrorFlag)

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
