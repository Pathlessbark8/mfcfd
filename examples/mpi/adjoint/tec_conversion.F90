program main
      implicit none

      integer::max_points,flag,i1,i2,i3,j,k,l,ele,i,count
      real*8::pr,rho,u1,u2,sen,x,y,eps1, eps2, mach, entropy
!
      real*8 :: xpc,xnc,ypc,ync
!
      integer :: flag1,qtdepth, max_depth, max_depth_point, count1
      real*8 :: sum, sd
      sum = 0.0d0
!

      open(101,file='output.dat')
      open(102,file='preprocessorfile.1.ele')
      open(103,file='meshfree.tec')
      open(104,file='sensor_flag.dat')
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
      write(103,*) 'VARIABLES = "X","Y","flag","qtdepth","rho","u1","u2","pr","mach","entropy","sensor"'
      write(103,*) 'ZONE T="P_1", DATAPACKING=POINT, NODES=',max_points,', ELEMENTS=',ele,', ZONETYPE=FETRIANGLE'

      do i=1,max_points
        read(101,*)x,y,flag,qtdepth,rho,u1,u2,pr,mach,entropy,sen
!        read(105,*)flag1,x,y,xpc,xnc,ypc,ync
        if(max_depth < qtdepth) then
                max_depth = qtdepth
                max_depth_point = i
        endif
!
        if(qtdepth .ge. 24) then
                print*, i, qtdepth
        endif
!        write(103,*)x,y,flag,qtdepth,rho,u1,u2,pr,mach,entropy,sen,xpc,xnc,ypc,ync
        write(103,*)x,y,flag,qtdepth,rho,u1,u2,pr,mach,entropy,sen
!        if(flag .ne. 0 .and. sen>eps1 .and. sen < eps2 ) then
        if(sen>eps1 .and. sen < eps2 .and. qtdepth .le. 18) then
                write(104,*)i,1
                count=count+1
        else if(sen < 1e-5 .and. flag .ne. 0 .and. flag .ne. 2) then
                write(104,*)i,2
                count1 = count1 + 1
        else
                        write(104,*)i,0
        endif
      sum = sum + sen*sen
      end do
!
!
      sd = dsqrt(sum/max_points)
!
!
      print*,'!!!!!!!!!!!!!!!!!!!!!'
      print*,'flag number:',count, count1, sd
      print*,'!!!!!!!!!!!!!!!!!!!!!'

      do i=1,ele
        read(102,*)j,i1,i2,i3
        write(103,*)i1,i2,i3
      end do
      

end program main
