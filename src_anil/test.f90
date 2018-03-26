subroutine test(h,x,y,z,m,n)
      implicit none
      integer :: h
      character(len=10) :: n 
      real*8 :: x,y,z,m
      
      open(h,file=n)

      write(h,'(5e30.20)') x,y,z,m


end subroutine



