!This subroutine is used for testing purposes. It will print out 
subroutine test(x,y,z,m,xx,yy,zz,mm)
        use data_structure_mod
      implicit none
      integer :: h
      real*8 :: x,y,z,m,xx,yy,zz,mm
      character(len=10) :: fi,itos

      h=rank+1

      fi=  'test/new'//trim(itos(1,rank))
      if(proc==1) fi=  'new'
      open(h,file=fi)

      write(h,'(8e30.20)') x,y,z,m,xx,yy,zz,mm


end subroutine



