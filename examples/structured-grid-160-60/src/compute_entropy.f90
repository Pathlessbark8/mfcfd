module compute_entropy_mod
!
	use data_structure_mod
!
	contains
!
!
		subroutine compute_entropy()
		!
		!
			implicit none
			
			integer :: k
			real*8 :: temp1, temp2
!			
!
			total_entropy = 0.d0
!			
			temp2 = dlog(pr_inf)
!
			do k = 1, max_points
					temp1 = (point(k)%rho)**gamma
					temp1 = point(k)%pr/temp1
					temp1 = dlog(temp1)
!				
					point(k)%entropy = dabs(temp1 - temp2)				
					total_entropy = total_entropy + dabs(temp1 - temp2)				
			enddo 
!
!			
		end subroutine !
!
!
end module compute_entropy_mod
!		
!
!
