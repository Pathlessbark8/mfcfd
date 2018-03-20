subroutine test(h,x,y,z)
        use data_structure_mod
      implicit none
      integer :: h
      real*8 :: x,y,z
      character(len=10) :: fi,itos

      fi=  'new'//trim(itos(1,rank))
      
      open(h,file=fi)

      write(h,'(5e30.20)') x,y,z


end subroutine



