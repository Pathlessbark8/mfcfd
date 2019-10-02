program main
      implicit none

      integer::max_points,flag,i1,i2,i3,j,k,l,ele,i,count
      real*8::pr,rho,u1,u2,sen,x,y,eps1, eps2, mach, entropy
!
!
      integer :: flag1,qtdepth, max_depth, max_depth_point, count1
      integer :: nx=160, ny=60
      real*8 :: sum, sd
      sum = 0.0d0
!

      open(101,file='output.dat')
      open(103,file='meshfree.dat')

      read(101,*)max_points

      write(103,*) 'TITLE="NACA 0012 case"'
      write(103,*) 'VARIABLES = "X","Y","flag","qtdepth","rho","u1","u2","pr","mach","entropy","sensor"'
      write(103,*) 'Zone I=',     nx, 'J=    ',ny, 'F=POINT'

      do i=1,max_points
        read(101,*)x,y,flag,qtdepth,rho,u1,u2,pr,mach,entropy,sen
        write(103,*)x,y,flag,qtdepth,rho,u1,u2,pr,mach,entropy,sen
      end do


end program main
