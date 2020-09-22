module point_preprocessor_mod

        use data_structure_mod_diff
        use petsc_data_structure_mod
        USE HDF5        

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
        allocate(point%conn(max_points,20))
        allocate(point%min_dist(max_points))
        allocate(point%left(max_points))
        allocate(point%right(max_points))
        allocate(point%original_id(local_points))

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

        do k = 1, local_points
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

            point%original_id(k) = H51DIntegerdataset(1)
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

            if(point%flag_2(k) > 0) then
                if(point%flag_2(k) > shapes)then
                    SETERRA(PETSC_COMM_WORLD,1,'shapes value wrong, check again')
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

        if (proc > 1) then
            allocate(pghost(ghost_points))
            ALLOCATE(H51DDoubledataset(3))

            ghost_group = trim(main_group)//'/ghost'
            CALL h5gopen_f(group_id, ghost_group, ghostgroup_id,ErrorFlag)
            IF (ErrorFlag.lt.0) THEN
                ErrorMessage=" *** Error opening ghost group "
                return
            ENDIF

            do k = 1, ghost_points
                point_string = '/'//trim(itos_unpad(k))
                dataset_string = trim(ghost_group)//point_string
    
                CALL h5dopen_f(ghostgroup_id , dataset_string , d_id, ErrorFlag)
                IF (ErrorFlag.lt.0) THEN
                    ErrorMessage=" *** Error opening dataset "
                    return
                ENDIF
                dims=1
                CALL h5dread_f(d_id, H5T_NATIVE_INTEGER, H5dataset, dims, ErrorFlag)
                pghost(k) = H5dataset
                CALL h5aopen_name_f(d_id, val1_s, a_id, ErrorFlag)
                IF (ErrorFlag.lt.0) THEN
                    ErrorMessage=" *** Error opening attribute "
                    return
                ENDIF
                dims=3
                CALL h5aread_f(a_id, H5T_NATIVE_DOUBLE, H51DDoubledataset, dims, ErrorFlag)
                CALL h5aclose_f(a_id,ErrorFlag)

                point%x(local_points + k) = H51DDoubledataset(1)
                point%y(local_points + k) = H51DDoubledataset(2)
                point%min_dist(local_points + k) = H51DDoubledataset(3)

                CALL h5dclose_f(d_id,ErrorFlag)
            end do

            deallocate(H51DDoubledataset)
            CALL h5gclose_f(ghostgroup_id, ErrorFlag)
            IF (ErrorFlag.lt.0) THEN
                ErrorMessage=" *** Error closing ghost group"
                return
            ENDIF
        end if

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


        ! subroutine read_phi_data()

        !         implicit none

        !         integer:: i, k, r, nproc, extra
        !         character(len=64) :: part_grid
        !         character(len=10) :: itos

        !         part_grid = 'phi_values_split/phi.dat'
        !         if (proc>1) part_grid = 'phi_values_split/phi-'//trim(itos(4,rank))//'.dat'

        !         OPEN(UNIT=101,FILE=trim(part_grid),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

        !         allocate(point%phi1(4, max_points))
        !         allocate(point%phi2(4, max_points))

        !         if (proc>1) then

        !                 do k = 1, local_points
                        
        !                         read(101,*) extra, point%phi1(1,k), point%phi1(2,k), point%phi1(3,k), &
        !                         & point%phi1(4,k) , point%phi2(1,k), point%phi2(2,k), point%phi2(3,k), &
        !                         & point%phi2(4,k)

        !                 end do

        !         else
        !                 do k = 1, max_points
                        
        !                         read(101,*) point%phi1(1,k), point%phi1(2,k), point%phi1(3,k), &
        !                         & point%phi1(4,k) , point%phi2(1,k), point%phi2(2,k), point%phi2(3,k), &
        !                         & point%phi2(4,k)

        !                 end do

        !         end if
                                        
        !         CLOSE(UNIT=101)
        
        ! end subroutine 

end module 
