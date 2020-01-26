module limiters_mod


        use data_structure_mod

contains


!	The following two subroutines are used for min-max limiter ..


        subroutine max_q_value(i, maxi)
        
                implicit none
                integer i, j, k, r
		real*8 :: maxi(4)
                
                maxi = point%q(:,i)
                
                do j = 1, point%nbhs(i)
                
                        k = point%conn(i,j)
        
                        do r = 1, 4
                                if( maxi(r) < point%q(r,k) ) then
                                        maxi(r) = point%q(r,k)
                                endif
                        enddo
                enddo
        
        end subroutine 



        subroutine min_q_value(i, mini)

                implicit none

                integer i, j, k, r
	        real*8 :: mini(4)
                
                mini = point%q(:,i)
        
                do j = 1, point%nbhs(i)
                
                        k = point%conn(i,j)

                        do r = 1, 4
                                if( mini(r) > point%q(r,k) ) then
                                        mini(r) = point%q(r,k)
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

                        q = point%q(r,k) 
                        del_neg = qtilde(r) - q
                        if(dabs(del_neg) .le. 10e-6) then
                                phi(r)=1.d0
       
                        else if(dabs(del_neg) .gt. 10e-6) then                   
                                if(del_neg .gt. 0.d0) then 
                                        del_pos = point%qm(1,r,k) - q
                                
                                else if(del_neg .lt. 0.d0) then
                                        del_pos = point%qm(2,r,k) - q
                                endif

                                epsi = VL_CONST*point%min_dist(k)
                                epsi = epsi**3.0d0

                                num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
                                num = num*del_neg + 2.0d0*del_neg*del_neg*del_pos

                                den = del_pos*del_pos + 2.0d0*del_neg*del_neg ! Denominator ..
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

end module

