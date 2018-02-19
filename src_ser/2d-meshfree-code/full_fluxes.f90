module full_fluxes_mod
contains
!
	subroutine flux_Gx(Gx, nx, ny, u1, u2, rho, pr)
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
		u1 = ut
		u2 = un
!
!
!     Expressions for the split fluxes ..	
!
		Gx(1) = rho*u1
!
		Gx(2) = pr + rho*u1*u1
!
		Gx(3) = rho*u1*u2
!
		temp1 = 0.5*(u1*u1 + u2*u2);
		rho_e = 2.5*pr + rho*temp1;
		Gx(4) = (pr + rho_e)*u1
! 
		end subroutine
!
!
!
	subroutine flux_Gy(Gy, nx, ny, u1, u2, rho, pr)
!
!
		implicit none
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
      	u1 = ut
      	u2 = un
!      
!     Expressions for the split fluxes ..	
!
      	Gy(1) = rho*u2
!
      	Gy(2) = rho*u1*u2
!
      	Gy(3) = pr + rho*u2*u2
!
      	temp1 = 0.5*(u1*u1 + u2*u2);
      	rho_e = 2.5*pr + rho*temp1;
      	Gy(4) = (pr + rho_e)*u2
! 
      end subroutine
!
!
!
end module full_fluxes_mod
