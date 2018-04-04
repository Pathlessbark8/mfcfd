module limiters_mod


        use data_structure_mod

contains


!	The following two subroutines are used for min-max limiter ..


        subroutine max_q_value(i, maxi)
        
                implicit none
                integer i, j, k, r
		real*8 :: maxi(4)
                
                maxi = p%q(:,i)
                
                do j = 1, p%nbhs(i)
                
                        k = p%conn(i,j)
        
                        do r = 1, 4
                                if( maxi(r) < p%q(r,k) ) then
                                        maxi(r) = p%q(r,k)
                                endif
                        enddo
                enddo
        
        end subroutine 



        subroutine min_q_value(i, mini)

                implicit none

                integer i, j, k, r
	        real*8 :: mini(4)
                
                mini = p%q(:,i)
        
                do j = 1, p%nbhs(i)
                
                        k = p%conn(i,j)

                        do r = 1, 4
                                if( mini(r) > p%q(r,k) ) then
                                        mini(r) = p%q(r,k)
                                endif
                        enddo
                enddo


        end subroutine 



!	The following subroutines are used for venkatakrishnan limiter .. 	



        subroutine venkat_limiter(qtilde, phi, k)


                implicit none
        
                integer :: r, k
		real*8 :: qtilde(4), phi(4)
		real*8 :: q, del_neg, del_pos
		real*8 :: max_q, min_q, ds, epsi, num, den, temp 

                do r = 1, 4

                        q = p%q(r,k) 
                        del_neg = qtilde(r) - q
                        if(dabs(del_neg) .le. 10e-6) then
                                phi(r)=1.d0
       
                        else if(dabs(del_neg) .gt. 10e-6) then                      
                                if(del_neg .gt. 0.d0) then 
                                        call maximum(k, r, max_q)
                                        del_pos = max_q - q
                                else if(del_neg .lt. 0.d0) then
                                        call minimum(k, r, min_q)               
                                        del_pos = min_q - q
                                endif

                                call smallest_dist(k, ds)

                                epsi = VL_CONST*ds
                                epsi = epsi**3.0

                                num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
                                num = num*del_neg + 2.0*del_neg*del_neg*del_pos

                                den = del_pos*del_pos + 2.0*del_neg*del_neg ! Denominator ..
                                den = den + del_neg*del_pos + epsi*epsi
                                den = den*del_neg

                                temp = num/den

                                if(temp .lt. 1.d0) then
                                        phi(r) = temp
                                else 
                                        phi(r) = 1.d0
                                endif

                        endif

                enddo

        end subroutine 



        subroutine maximum(k, r, max)

                implicit none 
                
                integer :: k, r, j, nbh
		real*8 :: max

                max = p%q(r,k)

                do j = 1, p%nbhs(k)
                        nbh = p%conn(k,j)

                        if(p%q(r,nbh) .gt. max) then 
                                max = p%q(r,nbh)
                        endif
                enddo
        end subroutine


        subroutine minimum(k, r, min)

                implicit none 

                integer :: k, r, j, nbh
		real*8 :: min

                min = p%q(r,k)

                do j = 1, p%nbhs(k)

                        nbh = p%conn(k,j)

                        if(p%q(r,nbh) .lt. min) then 
                                min = p%q(r,nbh)
                        endif
                enddo
        end subroutine!



        subroutine smallest_dist(k, min_dist)

                implicit none

                integer :: k, j, nbh
                real*8 :: dx, dy, ds, min_dist

                min_dist = 10000.d0

                do j = 1, p%nbhs(k)
                        nbh = p%conn(k,j)
                        dx = p%x(nbh) - p%x(k)
                        dy = p%y(nbh) - p%y(k)


                        ds = dsqrt(dx*dx + dy*dy)


                        if(ds .lt. min_dist) then
                                min_dist = ds
                         endif

                enddo   

   
        end subroutine!



end module

