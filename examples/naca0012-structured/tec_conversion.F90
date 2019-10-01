program main
      implicit none

      integer::max_points,flag,i1,i2,i3,j,k,l,ele,i,count
      real*8::pr,rho,u1,u2,sen,x,y,eps1, eps2, mach, entropy
!
!
      integer :: flag1,qtdepth, max_depth, max_depth_point, count1, dummy
      integer :: nx=160, ny=60
      integer :: parallel = 1
      real*8 :: objx,objy

      if(parallel == 0)open(101,file='solution/sol.dat')
      if(parallel == 1)open(101,file='sol.dat')
      open(103,file='meshfree.dat')

      if (parallel == 0) then

        write(103,*) 'TITLE="NACA 0012 case"'
        write(103,*) 'VARIABLES = "X","Y","objectivex","objectivey"'
        write(103,*) 'Zone I=',     nx, 'J=    ',ny, 'F=POINT'

        do i=1,nx*ny
                read(101,*)dummy,x,y,objx,objy
                write(103,*)x,y,objx,objy
        end do
      else
        call execute_command_line('sh merge.sh')
        write(103,*) 'TITLE="NACA 0012 case"'
        write(103,*) 'VARIABLES = "X","Y","objectivex","objectivey"'
        write(103,*) 'Zone I=',     nx, 'J=    ',ny, 'F=POINT'

        do i=1,nx*ny
                read(101,*)dummy,x,y,objx,objy
                write(103,*)x,y,objx,objy
        end do
              
      end if


end program main
