program main
      implicit none

      integer::max_points,flag,i1,i2,i3,j,k,l,ele,i,count
      real*8::pr,rho,u1,u2,sen,x,y,eps1, eps2, mach, entropy
!
!
      integer :: flag1,qtdepth, max_depth, max_depth_point, count1, dummy
      integer :: nx=160, ny=60
      integer :: parallel = 1
      real*8 :: sum, sd
      sum = 0.0d0
!

      if(parallel == 0)open(101,file='solution/sol.dat')
      if(parallel == 1)open(101,file='sol.dat')
      open(103,file='meshfree.dat')

      if (parallel == 0) then

        read(101,*)max_points

        write(103,*) 'TITLE="NACA 0012 case"'
        write(103,*) 'VARIABLES = "X","Y","rho","u1","u2","pr"'
        write(103,*) 'Zone I=',     nx, 'J=    ',ny, 'F=POINT'

        do i=1,max_points
                read(101,*)dummy,dummy,dummy,x,y,rho,u1,u2,pr
                write(103,*)x,y,rho,u1,u2,pr
        end do
      else
        call execute_command_line('sh merge.sh')
        write(103,*) 'TITLE="NACA 0012 case"'
        write(103,*) 'VARIABLES = "X","Y","rho","u1","u2","pr"'
        write(103,*) 'Zone I=',     nx, 'J=    ',ny, 'F=POINT'

        do i=1,nx*ny
                read(101,*)dummy,dummy,dummy,x,y,rho,u1,u2,pr
                write(103,*)x,y,rho,u1,u2,pr
        end do
              
      end if


end program main
