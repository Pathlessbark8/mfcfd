module miscellaneous_fluxes_mod
!
!
!
	contains
!
	
		SUBROUTINE flux_Gx(Gx, nx, ny, u1, u2, rho, pr)
!
!
			IMPLICIT NONE
!
!
			double precision Gx(4), u1, u2, rho, pr
			double precision tx, ty, nx, ny, ut, un
			double precision temp1, rho_e
! 
!
			tx = ny
			ty = -nx
!
			ut = u1*tx + u2*ty
			un = u1*nx + u2*ny
!
!      
!     Expressions for the split fluxes ..	
!
			Gx(1) = rho*ut
!
			Gx(2) = pr + rho*ut*ut
!
			Gx(3) = rho*ut*un
!
			temp1 = 0.5*(ut*ut + un*un);
			rho_e = 2.5*pr + rho*temp1;
			Gx(4) = (pr + rho_e)*ut
! 
	      end
!
!
!
!
      SUBROUTINE flux_Gy(Gy, nx, ny, u1, u2, rho, pr)
!
!
 	     IMPLICIT NONE
!
!
    	  double precision Gy(4), u1, u2, rho, pr
	      double precision tx, ty, nx, ny, ut, un
	      double precision temp1, rho_e
! 
!
	      tx = ny
	      ty = -nx
!
	      ut = u1*tx + u2*ty
	      un = u1*nx + u2*ny
!
!      
!     Expressions for the split fluxes ..	
!
	      Gy(1) = rho*un
!
	      Gy(2) = rho*ut*un
!
	      Gy(3) = pr + rho*un*un
!
	      temp1 = 0.5*(ut*ut + un*un);
	      rho_e = 2.5*pr + rho*temp1;
	      Gy(4) = (pr + rho_e)*un
! 
      end
!
!
!
end module miscellaneous_fluxes_mod