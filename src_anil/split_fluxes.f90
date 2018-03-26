module split_fluxes_mod
!
!
	use parameter_mod
!
contains
!
!
	subroutine flux_Gxp(Gxp, nx, ny, u1, u2, rho, pr)
!
!
		implicit none
!
		double precision Gxp(4), u1, u2, rho, pr
		double precision tx, ty, nx, ny, ut, un
		double precision beta, S1, B1, A1pos
		double precision temp1, temp2
		double precision pr_by_rho, u_sqr
! 
!
		tx = ny
		ty = -nx
!
		ut = u1*tx + u2*ty
		un = u1*nx + u2*ny
!
		beta = 0.5*rho/pr
		S1 = ut*dsqrt(beta) 
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		A1pos = 0.5*(1 + derf(S1))     
!
		pr_by_rho = pr/rho
		u_sqr = ut*ut + un*un
!
!     Expressions for the split fluxes ..	
!
		Gxp(1) = rho*(ut*A1pos + B1)  
!      
		temp1 = pr_by_rho + ut*ut
		temp2 = temp1*A1pos + ut*B1
		Gxp(2) = rho*temp2
!
		temp1 = ut*un*A1pos + un*B1
		Gxp(3) = rho*temp1
!
		temp1 = (7.0d0*pr_by_rho) + u_sqr
		temp2 = 0.5*ut*temp1*A1pos 
		temp1 = (6.0d0*pr_by_rho) + u_sqr
		Gxp(4) = rho*(temp2 + 0.5*temp1*B1)
!
!
      end
!
!
	subroutine flux_Gxn(Gxn, nx, ny, u1, u2, rho, pr)
!
!
		implicit none
!
		double precision Gxn(4), u1, u2, rho, pr
		double precision tx, ty, nx, ny, ut, un
		double precision beta, S1, B1, A1neg
		double precision temp1, temp2
		double precision pr_by_rho, u_sqr
! 
!
		tx = ny
		ty = -nx
!
		ut = u1*tx + u2*ty
		un = u1*nx + u2*ny
!
		beta = 0.5*rho/pr
		S1 = ut*dsqrt(beta) 
		B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
		A1neg = 0.5*(1 - derf(S1))     
!
		pr_by_rho = pr/rho
		u_sqr = ut*ut + un*un
!
!		Expressions for the split fluxes ..	
!
		Gxn(1) = rho*(ut*A1neg - B1)  
!
		temp1 = pr_by_rho + ut*ut
		temp2 = temp1*A1neg - ut*B1
		Gxn(2) = rho*temp2
!
		temp1 = ut*un*A1neg - un*B1
		Gxn(3) = rho*temp1
!
		temp1 = (7.0d0*pr_by_rho) + u_sqr
		temp2 = 0.5*ut*temp1*A1neg 
		temp1 = (6.0d0*pr_by_rho) + u_sqr
		Gxn(4) = rho*(temp2 - 0.5*temp1*B1)
!
!
      end	
!
!
	subroutine flux_Gyp(Gyp, nx, ny, u1, u2, rho, pr)
!
!
		implicit none
!
		double precision Gyp(4), u1, u2, rho, pr
		double precision tx, ty, nx, ny, ut, un
		double precision beta, S2, B2, A2pos
		double precision temp1, temp2
		double precision pr_by_rho, u_sqr
!
!
		tx = ny
		ty = -nx
!
		ut = u1*tx + u2*ty
		un = u1*nx + u2*ny
!
		beta = 0.5*rho/pr
		S2 = un*dsqrt(beta) 
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		A2pos = 0.5*(1 + derf(S2))     
!
		pr_by_rho = pr/rho
		u_sqr = ut*ut + un*un
!
!		Expressions for the split fluxes ..	
!
		Gyp(1) = rho*(un*A2pos + B2)  
!
		temp1 = pr_by_rho + un*un
		temp2 = temp1*A2pos + un*B2
		Gyp(3) = rho*temp2
!
		temp1 = ut*un*A2pos + ut*B2
		Gyp(2) = rho*temp1
!
		temp1 = (7.0d0*pr_by_rho) + u_sqr
		temp2 = 0.5*un*temp1*A2pos 
		temp1 = (6.0d0*pr_by_rho) + u_sqr
		Gyp(4) = rho*(temp2 + 0.5*temp1*B2)
!
!
      end
!
!
!
	subroutine flux_Gyn(Gyn, nx, ny, u1, u2, rho, pr)
!
!
		implicit none
!
!
		double precision Gyn(4), u1, u2, rho, pr
		double precision tx, ty, nx, ny, ut, un
		double precision beta, S2, B2, A2neg
		double precision temp1, temp2
		double precision pr_by_rho, u_sqr
!
!
		tx = ny
		ty = -nx
!
		ut = u1*tx + u2*ty
		un = u1*nx + u2*ny
!
		beta = 0.5*rho/pr
		S2 = un*dsqrt(beta) 
		B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
		A2neg = 0.5*(1 - derf(S2))     
!
		pr_by_rho = pr/rho
		u_sqr = ut*ut + un*un
!
!		Expressions for the split fluxes ..	
!
		Gyn(1) = rho*(un*A2neg - B2)  
! 
		temp1 = pr_by_rho + un*un
		temp2 = temp1*A2neg - un*B2
		Gyn(3) = rho*temp2
!
		temp1 = ut*un*A2neg - ut*B2
		Gyn(2) = rho*temp1
!
		temp1 = (7.0d0*pr_by_rho) + u_sqr
		temp2 = 0.5*un*temp1*A2neg 
		temp1 = (6.0d0*pr_by_rho) + u_sqr
		Gyn(4) = rho*(temp2 - 0.5*temp1*B2)
!
!
      end
!
!
end module split_fluxes_mod
