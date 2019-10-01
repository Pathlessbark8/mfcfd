program main
      implicit none

      integer::max_points,flag,i1,i2,i3,j,k,l,ele,i,count
      real*8::pr,rho,u1,u2,sen,x,y,eps1, eps2, mach, entropy, xb,yb
!
      real*8 :: xpc,xnc,ypc,ync, dummy
!
      integer :: flag1,qtdepth, max_depth, max_depth_point, count1
      real*8 :: sum, sd
      sum = 0.0d0
!

      open(101,file='sort.dat')
      open(102,file='preprocessorfile.1.ele')
      open(103,file='meshfree.tec')
!      open(104,file='sensor_flag.dat')
!      open(105,file='condition_numbers.dat')

      read(101,*)max_points
      read(102,*)ele,j,k

      eps1=0.10
      eps2 = 1.0
      count=0
      count1 = 0
      max_depth = 0

      write(103,*) 'TITLE="Subsonic flow over Williams airfoil"'
!      write(103,*) 'VARIABLES = "X","Y","flag","qtdepth","rho","u1","u2","pr","mach","entropy","sensor","xpc","xnc","ypc","ypn"'
      write(103,*) 'VARIABLES = "X","Y","xb","yb"'
      write(103,*) 'ZONE T="P_1", DATAPACKING=POINT, NODES=',max_points,', ELEMENTS=',ele,', ZONETYPE=FETRIANGLE'

      do i=1,max_points
        read(101,*)dummy,x,y,xb,yb
        write(103,*)x,y,xb,yb
!        if(flag .ne. 0 .and. sen>eps1 .and. sen < eps2 ) then
      end do
!
!
      do i=1,ele
        read(102,*)j,i1,i2,i3
        write(103,*)i1,i2,i3
      end do
      

end program main
