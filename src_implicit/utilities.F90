function itos(ndigit, n) result(str)
   implicit none
   integer,intent(in) :: ndigit, n
   character(len=10) :: str
   integer :: i, nmax, nd

   nmax = 10**ndigit - 1

   ! Count no. of zeroes to pad
   if(n <= 9) then
      write(str,'(I1)') n
      nd = ndigit - 1
   else if(n <= 99) then
      write(str,'(I2)') n
      nd = ndigit - 2
   else if(n <= 999) then
      write(str,'(I3)') n
      nd = ndigit - 3
   else if(n <= 9999) then
      write(str,'(I4)') n
      nd = ndigit - 4
   else
      nd = 0 ! Dummy to suppress compiler warning
      print*,'not recoginized'
   endif

   ! Pad with zeroes in beginning
   do i=1,nd
      str = '0'//trim(str)
   enddo

end function itos
