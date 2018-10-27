module post_processing_mod
	

        use data_structure_mod
	

contains


	subroutine print_primal_output()
	
	
		implicit none
			
		integer :: k,count,i
                count=0
		
		OPEN(UNIT=112,FILE='output.dat',FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		OPEN(UNIT=111,FILE="primal.poly",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		OPEN(UNIT=113,FILE="solution.dat",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
		
		!for traingle software
                write(111,*) local_points, 2 , 0 , 1
                do k=1, wall_points
!                        if(point%flag_1(k)==0) then
!                               write(112,*) point%x(k),','&4yy
!                               &,point%y(k),',',0,',',point%prim(1,k),',',point%prim(2,k),',',&
!                        &point%prim(3,k),',',point%prim(4,k)!,point%flux_res(1,k)&
! !                       &,point%flux_res(2,k),point%flux_res(3,k),point%flux_res(4,k)
!                        else
!
!                                write(111,*) point%x(k),','&
!                                &,point%y(k),',',0,',',point%prim(1,k),',',point%prim(2,k),',',&
!                        &point%prim(3,k),',',point%prim(4,k)!,point%flux_res(1,k)&
!                        endif
                                
                                write(111,*) k,point%x(k),point%y(k), 2
                enddo
                do k=wall_points+1,local_points
                          if(point%flag_1(k)==2) then
                                  write(111,*) k,point%x(k),point%y(k), 1
                                  count= count+1
                          else 
                                  write(111,*) k,point%x(k),point%y(k), 0
                          end if
                                  
                end do


                           write(111,*) wall_points+count,1
                do k=1,wall_points-1

                           write(111,*) k, k ,k+1,2
                end do

                write(111,*)wall_points,wall_points,1,2

                do k =wall_points+1,count+wall_points
                        i=outer_points_index(k-wall_points)
                           write(111,*) k, i,point%left(i),1
                           
                enddo

		write(111,*) 1
		write(111,*) 1,0.3161, 0.005348
!tecplot
                do k=1,max_points
                        
		        write(112,*)point%x(k),point%y(k),point%prim(1,k)&
                                &,point%prim(2,k),point%prim(3,k),point%prim(4,k),point%sensor(k)
                end do

                !solution data
                do k=1,max_points
                               write(113,*) k,point%x(k),&
                               &point%y(k),point%prim(1,k),point%prim(2,k),&
                        &point%prim(3,k),point%prim(4,k),point%flux_res(1,k)&
                        &,point%flux_res(2,k),point%flux_res(3,k),point%flux_res(4,k)
                end do
					
		CLOSE(UNIT=111)
		CLOSE(UNIT=112)
			
	end subroutine		
        subroutine writevtk()
                implicit none
                ! Local variables
                integer,parameter :: fid = 50
                integer       :: i, j, tc
                character(64) :: sfile
                character(len=10) :: int2str
                integer,save  :: fcount = 0

                sfile = 'output.vtk'

                open(fid,file=trim(sfile))
                write(fid,'("# vtk DataFile Version 3.0")')
                write(fid,'("freeflow")')
                write(fid,'("ASCII")')

                write(fid,'("DATASET POLYDATA")')
!   write(fid,'("FIELD FieldData 2")')
!   write(fid,'("TIME 1 1 double")')
!   write(fid,'(e15.8)') t
!   write(fid,'("CYCLE 1 1 int")')
!   write(fid,'(i8)') it

                write(fid,'("POINTS",i16," float")') local_points
                do i=1,local_points
                 write(fid,*)point%x(i),point%y(i), 0
                enddo

!   ! Order vertices according to vtk, see page 9 of
!   ! http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
!   write(fid,'("CELLS",i16,i16)') g%nc, 5*g%ntet+6*g%npyr+7*g%npri+9*g%nhex
!   do i=1,g%nc
!      tc = g%tcell(i)
!      write(fid,*) nvc(tc), (g%vcell(v2m(j,tc),i)-1,j=1,nvc(tc))
!   enddo
!
!   write(fid,'("CELL_TYPES",i16)') g%nc
!   do i=1,g%nc
!      tc = g%tcell(i)
!      write(fid,*) vtkcelltype(tc)
!   enddo
!
                write(fid,'("POINT_DATA",i16)') local_points
                
                write(fid,'("SCALARS ux float 1")')
                write(fid,'("LOOKUP_TABLE default")')
                do i=1,local_points
                 write(fid,*)point%prim(2,i)
                enddo
!
!                write(fid,'("SCALARS uy float 1")')
!                write(fid,'("LOOKUP_TABLE default")')
!                do i=1,local_points
!                 write(fid,*)point%prim(3,i)
!                enddo
!
!                write(fid,'("SCALARS rho float 1")')
!                write(fid,'("LOOKUP_TABLE default")')
!                do i=1,local_points
!                 write(fid,*)point%prim(1,i)
!                enddo
!
!                write(fid,'("SCALARS pressure float 1")')
!                write(fid,'("LOOKUP_TABLE default")')
!                do i=1,local_points
!                 write(fid,*)point%prim(4,i)
!                enddo
!                
!
!!                write(fid,'("SCALARS sensor1 float 1")')
!!                write(fid,'("LOOKUP_TABLE default")')
!!                do i=1,local_points
!!                 write(fid,*)point%sensor1(i)
!!                enddo
!!
!!                write(fid,'("SCALARS sensor2 float 1")')
!!                write(fid,'("LOOKUP_TABLE default")')
!!                do i=1,local_points
!!                 write(fid,*)point%sensor2(i)
!!                enddo
!!
                close(fid)


        end subroutine writevtk



end module post_processing_mod
