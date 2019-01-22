module quadrant_fluxes_mod

!	This module consists of quadrant split fluxes 
!	with respect to the x-coordinate direction ..

       
        use parameter_mod

contains


        subroutine flux_quad_GxI(G, nx, ny, u1, u2, rho, pr)


                implicit none

                double precision G(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta
                double precision S1, B1, S2, B2
                double precision A1neg, A2neg
                double precision temp1, temp2, temp3, temp4
                double precision pr_by_rho, u_sqr


                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5d0*rho/pr
                S1 = ut*dsqrt(beta) 
                S2 = un*dsqrt(beta) 
                B1 = 0.5d0*dexp(-S1*S1)/dsqrt(pi*beta)
                B2 = 0.5d0*dexp(-S2*S2)/dsqrt(pi*beta)
                A1neg = 0.5d0*(1.0d0 - derf(S1))     
                A2neg = 0.5d0*(1.0d0 - derf(S2))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!	Expressions for the split fluxes ..	

                G(1) = rho*A2neg*(ut*A1neg - B1)  

                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1neg - ut*B1
                G(2) = rho*A2neg*temp2

                temp1 = ut*A1neg - B1
                temp2 = un*A2neg - B2
                G(3) = rho*temp1*temp2

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5d0*ut*temp1*A1neg
 
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                temp3 = 0.5d0*B1*temp1 

                temp1 = ut*A1neg - B1
                temp4 = 0.5d0*rho*un*B2*temp1
      
                G(4) = rho*A2neg*(temp2 - temp3) - temp4
 

        end subroutine



        subroutine flux_quad_GxII(G, nx, ny, u1, u2, rho, pr)


                IMPLICIT NONE

                double precision G(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta
                double precision S1, B1, S2, B2
                double precision A1pos, A2neg
                double precision temp1, temp2, temp3, temp4
                double precision pr_by_rho, u_sqr


                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny


                beta = 0.5d0*rho/pr
                S1 = ut*dsqrt(beta) 
                S2 = un*dsqrt(beta) 
                B1 = 0.5d0*dexp(-S1*S1)/dsqrt(pi*beta)
                B2 = 0.5d0*dexp(-S2*S2)/dsqrt(pi*beta)
                A1pos = 0.5d0*(1.d0 + derf(S1))     
                A2neg = 0.5d0*(1.d0 - derf(S2))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!	Expressions for the split fluxes ..	

                G(1) = rho*A2neg*(ut*A1pos + B1)  

                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1pos + ut*B1
                G(2) = rho*A2neg*temp2

                temp1 = ut*A1pos + B1
                temp2 = un*A2neg - B2
                G(3) = rho*temp1*temp2

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5d0*ut*temp1*A1pos

                temp1 = (6.0d0*pr_by_rho) + u_sqr
                temp3 = 0.5d0*B1*temp1 

                temp1 = ut*A1pos + B1
                temp4 = 0.5d0*rho*un*B2*temp1

                G(4) = rho*A2neg*(temp2 + temp3) - temp4


        end subroutine



        subroutine flux_quad_GxIII(G, nx, ny, u1, u2, rho, pr)


        implicit none

                double precision G(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta
                double precision S1, B1, S2, B2
                double precision A1pos, A2pos
                double precision temp1, temp2, temp3, temp4
                double precision pr_by_rho, u_sqr


                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5d0*rho/pr
                S1 = ut*dsqrt(beta) 
                S2 = un*dsqrt(beta) 
                B1 = 0.5d0*dexp(-S1*S1)/dsqrt(pi*beta)
                B2 = 0.5d0*dexp(-S2*S2)/dsqrt(pi*beta)
                A1pos = 0.5d0*(1.0d0 + derf(S1))     
                A2pos = 0.5d0*(1.0d0 + derf(S2))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!	Expressions for the split fluxes ..	

                G(1) = rho*A2pos*(ut*A1pos + B1)  

                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1pos + ut*B1
                G(2) = rho*A2pos*temp2

                temp1 = ut*A1pos + B1
                temp2 = un*A2pos + B2
                G(3) = rho*temp1*temp2

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5d0*ut*temp1*A1pos

                temp1 = (6.0d0*pr_by_rho) + u_sqr
                temp3 = 0.5d0*B1*temp1 

                temp1 = ut*A1pos + B1
                temp4 = 0.5d0*rho*un*B2*temp1

                G(4) = rho*A2pos*(temp2 + temp3) + temp4


        end subroutine



      SUBROUTINE flux_quad_GxIV(G, nx, ny, u1, u2, rho, pr)


      IMPLICIT NONE

                double precision G(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta
                double precision S1, B1, S2, B2
                double precision A1neg, A2pos
                double precision temp1, temp2, temp3, temp4
                double precision pr_by_rho, u_sqr
           
          
                tx = ny
                ty = -nx
          
                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny
          
                beta = 0.5d0*rho/pr
                S1 = ut*dsqrt(beta) 
                S2 = un*dsqrt(beta) 
                B1 = 0.5d0*dexp(-S1*S1)/dsqrt(pi*beta)
                B2 = 0.5d0*dexp(-S2*S2)/dsqrt(pi*beta)
                A1neg = 0.5d0*(1.0d0 - derf(S1))     
                A2pos = 0.5d0*(1.0d0 + derf(S2))     
          
                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un
          
          !	Expressions for the split fluxes ..	
          
                G(1) = rho*A2pos*(ut*A1neg - B1)  
                 
                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1neg - ut*B1
                G(2) = rho*A2pos*temp2
          
                temp1 = ut*A1neg - B1
                temp2 = un*A2pos + B2
                G(3) = rho*temp1*temp2
          
                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5d0*ut*temp1*A1neg
          
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                temp3 = 0.5d0*B1*temp1 
          
                temp1 = ut*A1neg - B1
                temp4 = 0.5d0*rho*un*B2*temp1
                
                G(4) = rho*A2pos*(temp2 - temp3) + temp4
           
        end

      !        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
!
!  Differentiation of flux_quad_gxi in forward (tangent) mode:
!   variations   of useful results: g
!   with respect to varying inputs: u
!   RW status of diff variables: g:out u:in
!
        SUBROUTINE FLUX_QUAD_GXI_D(g, gd, u, ud, nx, ny)
          IMPLICIT NONE
         
        
       
          DOUBLE PRECISION :: g(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta
          DOUBLE PRECISION :: betad
          DOUBLE PRECISION :: s1, b1, s2, b2
          DOUBLE PRECISION :: s1d, b1d, s2d, b2d
          DOUBLE PRECISION :: a1neg, a2neg
          DOUBLE PRECISION :: a1negd, a2negd
          DOUBLE PRECISION :: temp1, temp2, temp3, temp4
          DOUBLE PRECISION :: temp1d, temp2d, temp3d, temp4d
          DOUBLE PRECISION :: pr_by_rho, u_sqr
          DOUBLE PRECISION :: pr_by_rhod, u_sqrd
        !  INTRINSIC DSQRT
        !  INTRINSIC DEXP
        !  EXTERNAL DERF
        !  EXTERNAL DERF_D
        !  REAL :: DERF
        !  REAL :: DERF_D
          DOUBLE PRECISION :: result1
          DOUBLE PRECISION :: result1d
          DOUBLE PRECISION :: arg1
          DOUBLE PRECISION :: arg1d
          REAL*8 :: result10
          REAL*8 :: result10d
        ! 
        !
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5d0*(2.0d0*u(2)*ud(2)+2*u(3)*ud(3))/u(1)+0.5d0*ud(1)*(u(2)**2+u&
        &   (3)**2)/u(1)**2)/2.5d0
          pr = (u(4)-0.5d0/u(1)*(u(2)**2+u(3)**2))/2.5d0
        
        !
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
        !
          betad = (0.5d0*rhod*pr-0.5d0*rho*prd)/pr**2
          beta = 0.5d0*rho/pr
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5d0*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5d0*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5d0*DEXP(-(s1*s1))/result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5d0*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5d0*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5d0*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s1, s1d, result10)
        
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
        
          a1negd = -(0.5d0*result10d)
          a1neg = 0.5d0*(1.0d0-result10)
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2negd = -(0.5d0*result10d)
          a2neg = 0.5d0*(1.0d0-result10)
        !
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !
        !	Expressions for the split fluxes ..	
        !
          gd = 0.D0
          gd(1) = (rhod*a2neg+rho*a2negd)*(ut*a1neg-b1) + rho*a2neg*(utd*a1neg+&
        &   ut*a1negd-b1d)
          g(1) = rho*a2neg*(ut*a1neg-b1)
        !
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1neg + temp1*a1negd - utd*b1 - ut*b1d
          temp2 = temp1*a1neg - ut*b1
          gd(2) = (rhod*a2neg+rho*a2negd)*temp2 + rho*a2neg*temp2d
          g(2) = rho*a2neg*temp2
        !
          temp1d = utd*a1neg + ut*a1negd - b1d
          temp1 = ut*a1neg - b1
          temp2d = und*a2neg + un*a2negd - b2d
          temp2 = un*a2neg - b2
          gd(3) = (rhod*temp1+rho*temp1d)*temp2 + rho*temp1*temp2d
          g(3) = rho*temp1*temp2
        !
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5d0*((utd*temp1+ut*temp1d)*a1neg+ut*temp1*a1negd)
          temp2 = 0.5d0*ut*temp1*a1neg
        ! 
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          temp3d = 0.5d0*(b1d*temp1+b1*temp1d)
          temp3 = 0.5d0*b1*temp1
        !
          temp1d = utd*a1neg + ut*a1negd - b1d
          temp1 = ut*a1neg - b1
          temp4d = 0.5d0*((rhod*un+rho*und)*b2*temp1+rho*un*(b2d*temp1+b2*temp1d))
          temp4 = 0.5d0*rho*un*b2*temp1
        !      
          gd(4) = (rhod*a2neg+rho*a2negd)*(temp2-temp3) + rho*a2neg*(temp2d-&
        &   temp3d) - temp4d
          g(4) = rho*a2neg*(temp2-temp3) - temp4
        END SUBROUTINE FLUX_QUAD_GXI_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_quad_gxii in forward (tangent) mode:
        !   variations   of useful results: g
        !   with respect to varying inputs: u
        !   RW status of diff variables: g:out u:in
        !
        !
        !
        SUBROUTINE FLUX_QUAD_GXII_D(g, gd, u, ud, nx, ny)
          IMPLICIT NONE
        !
        !
        !
          DOUBLE PRECISION :: g(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta
          DOUBLE PRECISION :: betad
          DOUBLE PRECISION :: s1, b1, s2, b2
          DOUBLE PRECISION :: s1d, b1d, s2d, b2d
          DOUBLE PRECISION :: a1pos, a2neg
          DOUBLE PRECISION :: a1posd, a2negd
          DOUBLE PRECISION :: temp1, temp2, temp3, temp4
          DOUBLE PRECISION :: temp1d, temp2d, temp3d, temp4d
          DOUBLE PRECISION :: pr_by_rho, u_sqr
          DOUBLE PRECISION :: pr_by_rhod, u_sqrd
        !  INTRINSIC DSQRT
        !  INTRINSIC DEXP
        !  EXTERNAL DERF
        !  EXTERNAL DERF_D
        !  REAL :: DERF
        !  REAL :: DERF_D
          DOUBLE PRECISION :: result1
          DOUBLE PRECISION :: result1d
          DOUBLE PRECISION :: arg1
          DOUBLE PRECISION :: arg1d
          REAL*8 :: result10
          REAL*8 :: result10d
        ! 
        !
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5d0*(2.0d0*u(2)*ud(2)+2.0d0*u(3)*ud(3))/u(1)+0.5d0*ud(1)*(u(2)**2+u&
        &   (3)**2)/u(1)**2)/2.5d0
          pr = (u(4)-0.5d0/u(1)*(u(2)**2+u(3)**2))/2.5d0
          
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
        !
          betad = (0.5d0*rhod*pr-0.5d0*rho*prd)/pr**2
          beta = 0.5d0*rho/pr
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5d0*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5d0*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5d0*DEXP(-(s1*s1))/result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5d0*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5d0*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5d0*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s1, s1d, result10)
        
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
        
          a1posd = 0.5d0*result10d
          a1pos = 0.5d0*(1.0d0+result10)
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2negd = -(0.5d0*result10d)
          a2neg = 0.5d0*(1.0d0-result10)
        !
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !
        !	Expressions for the split fluxes ..	
        !
          gd = 0.D0
          gd(1) = (rhod*a2neg+rho*a2negd)*(ut*a1pos+b1) + rho*a2neg*(utd*a1pos+&
        &   ut*a1posd+b1d)
          g(1) = rho*a2neg*(ut*a1pos+b1)
        !
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1pos + temp1*a1posd + utd*b1 + ut*b1d
          temp2 = temp1*a1pos + ut*b1
          gd(2) = (rhod*a2neg+rho*a2negd)*temp2 + rho*a2neg*temp2d
          g(2) = rho*a2neg*temp2
        !
          temp1d = utd*a1pos + ut*a1posd + b1d
          temp1 = ut*a1pos + b1
          temp2d = und*a2neg + un*a2negd - b2d
          temp2 = un*a2neg - b2
          gd(3) = (rhod*temp1+rho*temp1d)*temp2 + rho*temp1*temp2d
          g(3) = rho*temp1*temp2
        !
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5d0*((utd*temp1+ut*temp1d)*a1pos+ut*temp1*a1posd)
          temp2 = 0.5d0*ut*temp1*a1pos
        !
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          temp3d = 0.5d0*(b1d*temp1+b1*temp1d)
          temp3 = 0.5d0*b1*temp1
        !
          temp1d = utd*a1pos + ut*a1posd + b1d
          temp1 = ut*a1pos + b1
          temp4d = 0.5d0*((rhod*un+rho*und)*b2*temp1+rho*un*(b2d*temp1+b2*temp1d))
          temp4 = 0.5d0*rho*un*b2*temp1
        !
          gd(4) = (rhod*a2neg+rho*a2negd)*(temp2+temp3) + rho*a2neg*(temp2d+&
        &   temp3d) - temp4d
          g(4) = rho*a2neg*(temp2+temp3) - temp4
        END SUBROUTINE FLUX_QUAD_GXII_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_quad_gxiii in forward (tangent) mode:
        !   variations   of useful results: g
        !   with respect to varying inputs: u
        !   RW status of diff variables: g:out u:in
        !
        !
        !
        SUBROUTINE FLUX_QUAD_GXIII_D(g, gd, u, ud, nx, ny)
          IMPLICIT NONE
        ! 
        !
        !
          DOUBLE PRECISION :: g(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta
          DOUBLE PRECISION :: betad
          DOUBLE PRECISION :: s1, b1, s2, b2
          DOUBLE PRECISION :: s1d, b1d, s2d, b2d
          DOUBLE PRECISION :: a1pos, a2pos
          DOUBLE PRECISION :: a1posd, a2posd
          DOUBLE PRECISION :: temp1, temp2, temp3, temp4
          DOUBLE PRECISION :: temp1d, temp2d, temp3d, temp4d
          DOUBLE PRECISION :: pr_by_rho, u_sqr
          DOUBLE PRECISION :: pr_by_rhod, u_sqrd
        !  INTRINSIC DSQRT
        !  INTRINSIC DEXP
        !  EXTERNAL DERF
        !  EXTERNAL DERF_D
        !  REAL :: DERF
        !  REAL :: DERF_D
          DOUBLE PRECISION :: result1
          DOUBLE PRECISION :: result1d
          DOUBLE PRECISION :: arg1
          DOUBLE PRECISION :: arg1d
          REAL*8 :: result10
          REAL*8 :: result10d
        ! 
        !
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5d0*(2.0d0*u(2)*ud(2)+2.0d0*u(3)*ud(3))/u(1)+0.5d0*ud(1)*(u(2)**2+u&
        &   (3)**2)/u(1)**2)/2.5d0
          pr = (u(4)-0.5d0/u(1)*(u(2)**2+u(3)**2))/2.5d0
        !
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
        !
          betad = (0.5d0*rhod*pr-0.5d0*rho*prd)/pr**2
          beta = 0.5d0*rho/pr
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5d0*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5d0*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5d0*DEXP(-(s1*s1))/result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5d0*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5d0*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5d0*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s1, s1d, result10)
          
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
        
          a1posd = 0.5d0*result10d
          a1pos = 0.5d0*(1.0d0+result10)
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2posd = 0.5d0*result10d
          a2pos = 0.5d0*(1.0d0+result10)
        !
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !
        !	Expressions for the split fluxes ..	
        !
          gd = 0.D0
          gd(1) = (rhod*a2pos+rho*a2posd)*(ut*a1pos+b1) + rho*a2pos*(utd*a1pos+&
        &   ut*a1posd+b1d)
          g(1) = rho*a2pos*(ut*a1pos+b1)
        !
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1pos + temp1*a1posd + utd*b1 + ut*b1d
          temp2 = temp1*a1pos + ut*b1
          gd(2) = (rhod*a2pos+rho*a2posd)*temp2 + rho*a2pos*temp2d
          g(2) = rho*a2pos*temp2
        !
          temp1d = utd*a1pos + ut*a1posd + b1d
          temp1 = ut*a1pos + b1
          temp2d = und*a2pos + un*a2posd + b2d
          temp2 = un*a2pos + b2
          gd(3) = (rhod*temp1+rho*temp1d)*temp2 + rho*temp1*temp2d
          g(3) = rho*temp1*temp2
        !
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5d0*((utd*temp1+ut*temp1d)*a1pos+ut*temp1*a1posd)
          temp2 = 0.5d0*ut*temp1*a1pos
        !
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          temp3d = 0.5d0*(b1d*temp1+b1*temp1d)
          temp3 = 0.5d0*b1*temp1
        !
          temp1d = utd*a1pos + ut*a1posd + b1d
          temp1 = ut*a1pos + b1
          temp4d = 0.5d0*((rhod*un+rho*und)*b2*temp1+rho*un*(b2d*temp1+b2*temp1d))
          temp4 = 0.5d0*rho*un*b2*temp1
        !
          gd(4) = (rhod*a2pos+rho*a2posd)*(temp2+temp3) + rho*a2pos*(temp2d+&
        &   temp3d) + temp4d
          g(4) = rho*a2pos*(temp2+temp3) + temp4
        END SUBROUTINE FLUX_QUAD_GXIII_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_quad_gxiv in forward (tangent) mode:
        !   variations   of useful results: g
        !   with respect to varying inputs: u
        !   RW status of diff variables: g:out u:in
        !
        !
        !
        SUBROUTINE FLUX_QUAD_GXIV_D(g, gd, u, ud, nx, ny)
          IMPLICIT NONE
        ! 
        !
          DOUBLE PRECISION :: g(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta
          DOUBLE PRECISION :: betad
          DOUBLE PRECISION :: s1, b1, s2, b2
          DOUBLE PRECISION :: s1d, b1d, s2d, b2d
          DOUBLE PRECISION :: a1neg, a2pos
          DOUBLE PRECISION :: a1negd, a2posd
          DOUBLE PRECISION :: temp1, temp2, temp3, temp4
          DOUBLE PRECISION :: temp1d, temp2d, temp3d, temp4d
          DOUBLE PRECISION :: pr_by_rho, u_sqr
          DOUBLE PRECISION :: pr_by_rhod, u_sqrd
        !  INTRINSIC DSQRT
        !  INTRINSIC DEXP
        !  EXTERNAL DERF
        !  EXTERNAL DERF_D
        !  REAL :: DERF
        !  REAL :: DERF_D
          DOUBLE PRECISION :: result1
          DOUBLE PRECISION :: result1d
          DOUBLE PRECISION :: arg1
          DOUBLE PRECISION :: arg1d
          REAL*8 :: result10
          REAL*8 :: result10d
        ! 
        !
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5d0*(2.0d0*u(2)*ud(2)+2.0d0*u(3)*ud(3))/u(1)+0.5d0*ud(1)*(u(2)**2+u&
        &   (3)**2)/u(1)**2)/2.5d0
          pr = (u(4)-0.5d0/u(1)*(u(2)**2+u(3)**2))/2.5d0
        !
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
        !
          betad = (0.5d0*rhod*pr-0.5d0*rho*prd)/pr**2
          beta = 0.5d0*rho/pr
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          IF (beta .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5d0*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5d0*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5d0*DEXP(-(s1*s1))/result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0d0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5d0*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5d0*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5d0*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s1, s1d, result10)
        
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
        
          a1negd = -(0.5d0*result10d)
          a1neg = 0.5d0*(1.0d0-result10)
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2posd = 0.5d0*result10d
          a2pos = 0.5d0*(1.0d0+result10)
        !
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !
        !	Expressions for the split fluxes ..	
        !
          gd = 0.D0
          gd(1) = (rhod*a2pos+rho*a2posd)*(ut*a1neg-b1) + rho*a2pos*(utd*a1neg+&
        &   ut*a1negd-b1d)
          g(1) = rho*a2pos*(ut*a1neg-b1)
        !       
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1neg + temp1*a1negd - utd*b1 - ut*b1d
          temp2 = temp1*a1neg - ut*b1
          gd(2) = (rhod*a2pos+rho*a2posd)*temp2 + rho*a2pos*temp2d
          g(2) = rho*a2pos*temp2
        !
          temp1d = utd*a1neg + ut*a1negd - b1d
          temp1 = ut*a1neg - b1
          temp2d = und*a2pos + un*a2posd + b2d
          temp2 = un*a2pos + b2
          gd(3) = (rhod*temp1+rho*temp1d)*temp2 + rho*temp1*temp2d
          g(3) = rho*temp1*temp2
        !
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5d0*((utd*temp1+ut*temp1d)*a1neg+ut*temp1*a1negd)
          temp2 = 0.5d0*ut*temp1*a1neg
        !
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          temp3d = 0.5d0*(b1d*temp1+b1*temp1d)
          temp3 = 0.5d0*b1*temp1
        !
          temp1d = utd*a1neg + ut*a1negd - b1d
          temp1 = ut*a1neg - b1
          temp4d = 0.5d0*((rhod*un+rho*und)*b2*temp1+rho*un*(b2d*temp1+b2*temp1d))
          temp4 = 0.5d0*rho*un*b2*temp1
        !      
          gd(4) = (rhod*a2pos+rho*a2posd)*(temp2-temp3) + rho*a2pos*(temp2d-&
        &   temp3d) + temp4d
          g(4) = rho*a2pos*(temp2-temp3) + temp4
        END SUBROUTINE FLUX_QUAD_GXIV_D

!
!
end module quadrant_fluxes_mod	
