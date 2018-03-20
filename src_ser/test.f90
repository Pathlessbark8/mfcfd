subroutine test(h,x,y,z,n)
      implicit none
      integer :: h
      character(len=10) :: n 
      real*8 :: x,y,z
      
      open(h,file=n)

      write(h,'(5e30.20)') x,y,z


end subroutine



