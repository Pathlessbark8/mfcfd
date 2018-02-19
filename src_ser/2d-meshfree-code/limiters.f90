module limiters_mod
!
!
	use data_structure_mod
!
contains
!
!
!	The following two subroutines are used for min-max limiter ..
!
!
	subroutine max_q_value(i, maxi)
	
		implicit none

		integer i, j, k, r
		real*8 :: maxi(4)
		
		maxi = point(i)%q
		
		do j = 1, point(i)%nbhs
		
				k = point(i)%conn(j)
				
				do r = 1, 4
					if( maxi(r) < point(k)%q(r) ) then
							maxi(r) = point(k)%q(r)
					endif
				enddo	
		enddo
		
				
	end subroutine 
!
!
!		
	subroutine min_q_value(i, mini)
	
		implicit none

		integer i, j, k, r
		real*8 :: mini(4)
		
		mini = point(i)%q
		
		do j = 1, point(i)%nbhs
		
				k = point(i)%conn(j)
				
				do r = 1, 4
					if( mini(r) > point(k)%q(r) ) then
							mini(r) = point(k)%q(r)
					endif
				enddo	
		enddo
		
		
	end subroutine 
!
!
!
!	The following subroutines are used for venkatakrishnan limiter .. 	
!
!
!
		subroutine venkat_limiter(qtilde, phi, k)
!
!
			implicit none
			
			integer :: r, k
			real*8 :: qtilde(4), phi(4)
			real*8 :: q, del_neg, del_pos
			real*8 :: max_q, min_q, ds, epsi, num, den, temp 
		
				do r = 1, 4
				
						q = point(k)%q(r) 
		                del_neg = qtilde(r) - q

		                if(dabs(del_neg) .le. 10e-6) then
		                       phi(r)=1.d0
!		                       
						else if(dabs(del_neg) .gt. 10e-6) then		                       
								if(del_neg .gt. 0.d0) then 
										call maximum(k, r, max_q)
										del_pos = max_q - q
		                        else if(del_neg .lt. 0.d0) then
										call minimum(k, r, min_q)		                        
	        	                        del_pos = min_q - q
								endif
!
								call smallest_dist(k, ds)
!
								epsi = VL_CONST*ds
		                        epsi = epsi**3.0
!
 								num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
								num = num*del_neg + 2.0*del_neg*del_neg*del_pos
!
 								den = del_pos*del_pos + 2.0*del_neg*del_neg ! Denominator ..
								den = den + del_neg*del_pos + epsi*epsi
								den = den*del_neg
!
								temp = num/den
!
								if(temp .lt. 1.d0) then
										phi(r) = temp
								else 
										phi(r) = 1.d0
								endif	
!
						endif								
!
				enddo
!   		
			end subroutine !			
!
!
!
			subroutine maximum(k, r, max)
			!
				implicit none 
				!
					integer :: k, r, j, nbh
					real*8 :: max
!	
					max = point(k)%q(r)
!
					do j = 1, point(k)%nbhs					
							nbh = point(k)%conn(j)
!
							if(point(nbh)%q(r) .gt. max) then 
									max = point(nbh)%q(r)
							endif
					enddo		
			end subroutine!
!
!
			subroutine minimum(k, r, min)
			!
				implicit none 
				!
					integer :: k, r, j, nbh
					real*8 :: min
!
					min = point(k)%q(r)
!
					do j = 1, point(k)%nbhs					
							nbh = point(k)%conn(j)
!
							if(point(nbh)%q(r) .lt. min) then 
									min = point(nbh)%q(r)
							endif
					enddo		
			end subroutine!
!
!
!
			subroutine smallest_dist(k, min_dist)
!
			implicit none	
				!
				integer :: k, j, nbh
				real*8 :: dx, dy, ds, min_dist
!
					min_dist = 10000.d0

					do j = 1, point(k)%nbhs					
							nbh = point(k)%conn(j)
							
							dx = point(nbh)%x - point(k)%x
							dy = point(nbh)%y - point(k)%y
!
!
 							ds = dsqrt(dx*dx + dy*dy)
!
!
				            if(ds .lt. min_dist) then
				            		min_dist = ds
				            endif
!
					enddo    		
!
!   
			end subroutine!
!
!
!		
end module
		
				
												
							
						
						
