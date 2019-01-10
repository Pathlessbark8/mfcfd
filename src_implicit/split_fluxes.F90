module split_fluxes_mod


        use parameter_mod

contains

        subroutine flux_Gxp(Gxp, nx, ny, u1, u2, rho, pr)

                implicit none

                double precision Gxp(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta, S1, B1, A1pos
                double precision temp1, temp2
                double precision pr_by_rho, u_sqr


                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5*rho/pr
                S1 = ut*dsqrt(beta) 
                B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
                A1pos = 0.5*(1 + derf(S1))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!     Expressions for the split fluxes ..	

                Gxp(1) = rho*(ut*A1pos + B1)  

                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1pos + ut*B1
                Gxp(2) = rho*temp2

                temp1 = ut*un*A1pos + un*B1
                Gxp(3) = rho*temp1

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5*ut*temp1*A1pos 
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                Gxp(4) = rho*(temp2 + 0.5*temp1*B1)


        end


        subroutine flux_Gxn(Gxn, nx, ny, u1, u2, rho, pr)


                implicit none

                double precision Gxn(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta, S1, B1, A1neg
                double precision temp1, temp2
                double precision pr_by_rho, u_sqr

                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5*rho/pr
                S1 = ut*dsqrt(beta) 
                B1 = 0.5*dexp(-S1*S1)/dsqrt(pi*beta)
                A1neg = 0.5*(1 - derf(S1))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!		Expressions for the split fluxes ..	

                Gxn(1) = rho*(ut*A1neg - B1)  

                temp1 = pr_by_rho + ut*ut
                temp2 = temp1*A1neg - ut*B1
                Gxn(2) = rho*temp2

                temp1 = ut*un*A1neg - un*B1
                Gxn(3) = rho*temp1

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5*ut*temp1*A1neg 
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                Gxn(4) = rho*(temp2 - 0.5*temp1*B1)


      end


        subroutine flux_Gyp(Gyp, nx, ny, u1, u2, rho, pr)


                implicit none

                double precision Gyp(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta, S2, B2, A2pos
                double precision temp1, temp2
                double precision pr_by_rho, u_sqr

                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5*rho/pr
                S2 = un*dsqrt(beta) 
                B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
                A2pos = 0.5*(1 + derf(S2))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!		Expressions for the split fluxes ..	

                Gyp(1) = rho*(un*A2pos + B2)  
        
                temp1 = pr_by_rho + un*un
                temp2 = temp1*A2pos + un*B2
                Gyp(3) = rho*temp2

                temp1 = ut*un*A2pos + ut*B2
                Gyp(2) = rho*temp1

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5*un*temp1*A2pos 
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                Gyp(4) = rho*(temp2 + 0.5*temp1*B2)


      end


        subroutine flux_Gyn(Gyn, nx, ny, u1, u2, rho, pr)


                implicit none


                double precision Gyn(4), u1, u2, rho, pr
                double precision tx, ty, nx, ny, ut, un
                double precision beta, S2, B2, A2neg
                double precision temp1, temp2
                double precision pr_by_rho, u_sqr


                tx = ny
                ty = -nx

                ut = u1*tx + u2*ty
                un = u1*nx + u2*ny

                beta = 0.5*rho/pr
                S2 = un*dsqrt(beta) 
                B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*beta)
                A2neg = 0.5*(1 - derf(S2))     

                pr_by_rho = pr/rho
                u_sqr = ut*ut + un*un

!		Expressions for the split fluxes ..	

                Gyn(1) = rho*(un*A2neg - B2)  

                temp1 = pr_by_rho + un*un
                temp2 = temp1*A2neg - un*B2
                Gyn(3) = rho*temp2

                temp1 = ut*un*A2neg - ut*B2
                Gyn(2) = rho*temp1

                temp1 = (7.0d0*pr_by_rho) + u_sqr
                temp2 = 0.5*un*temp1*A2neg 
                temp1 = (6.0d0*pr_by_rho) + u_sqr
                Gyn(4) = rho*(temp2 - 0.5*temp1*B2)


      end

              !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_gxp in forward (tangent) mode:
        !   variations   of useful results: gxp
        !   with respect to varying inputs: u
        !   RW status of diff variables: gxp:out u:in
        SUBROUTINE FLUX_GXP_D(gxp, gxpd, u, ud, nx, ny)
          IMPLICIT NONE
          DOUBLE PRECISION :: gxp(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gxpd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta, s1, b1, a1pos
          DOUBLE PRECISION :: betad, s1d, b1d, a1posd
          DOUBLE PRECISION :: temp1, temp2
          DOUBLE PRECISION :: temp1d, temp2d
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
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5*(2*u(2)*ud(2)+2*u(3)*ud(3))/u(1)+0.5*ud(1)*(u(2)**2+u&
&   (3)**2)/u(1)**2)/2.5
          pr = (u(4)-0.5/u(1)*(u(2)**2+u(3)**2))/2.5
!
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
          betad = (0.5*rhod*pr-0.5*rho*prd)/pr**2
          beta = 0.5*rho/pr
          IF (beta .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5*DEXP(-(s1*s1))/result1
        !  result10d = DERF_D(s1, s1d, result10)
          
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
          
          a1posd = 0.5*result10d
          a1pos = 0.5*(1+result10)
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !     Expressions for the split fluxes ..	
          gxpd = 0.D0
          gxpd(1) = rhod*(ut*a1pos+b1) + rho*(utd*a1pos+ut*a1posd+b1d)
          gxp(1) = rho*(ut*a1pos+b1)
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1pos + temp1*a1posd + utd*b1 + ut*b1d
          temp2 = temp1*a1pos + ut*b1
          gxpd(2) = rhod*temp2 + rho*temp2d
          gxp(2) = rho*temp2
          temp1d = (utd*un+ut*und)*a1pos + ut*un*a1posd + und*b1 + un*b1d
          temp1 = ut*un*a1pos + un*b1
          gxpd(3) = rhod*temp1 + rho*temp1d
          gxp(3) = rho*temp1
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5*((utd*temp1+ut*temp1d)*a1pos+ut*temp1*a1posd)
          temp2 = 0.5*ut*temp1*a1pos
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          gxpd(4) = rhod*(temp2+0.5*temp1*b1) + rho*(temp2d+0.5*(temp1d*b1+temp1&
        &   *b1d))
          gxp(4) = rho*(temp2+0.5*temp1*b1)
        END SUBROUTINE FLUX_GXP_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_gxn in forward (tangent) mode:
        !   variations   of useful results: gxn
        !   with respect to varying inputs: u
        !   RW status of diff variables: gxn:out u:in
        SUBROUTINE FLUX_GXN_D(gxn, gxnd, u, ud, nx, ny)
          IMPLICIT NONE
          DOUBLE PRECISION :: gxn(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gxnd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta, s1, b1, a1neg
          DOUBLE PRECISION :: betad, s1d, b1d, a1negd
          DOUBLE PRECISION :: temp1, temp2
          DOUBLE PRECISION :: temp1d, temp2d
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
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5*(2*u(2)*ud(2)+2*u(3)*ud(3))/u(1)+0.5*ud(1)*(u(2)**2+u&
&   (3)**2)/u(1)**2)/2.5
          pr = (u(4)-0.5/u(1)*(u(2)**2+u(3)**2))/2.5
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
          betad = (0.5*rhod*pr-0.5*rho*prd)/pr**2
          beta = 0.5*rho/pr
          IF (beta .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s1d = utd*result1 + ut*result1d
          s1 = ut*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b1d = (-(0.5*(s1d*s1+s1*s1d)*DEXP(-(s1*s1))*result1)-0.5*DEXP(-(s1*s1)&
        &   )*result1d)/result1**2
          b1 = 0.5*DEXP(-(s1*s1))/result1
        !  result10d = DERF_D(s1, s1d, result10)
          
          result10d = dexp(-s1**2)*(2.d0/sqrt(pi))*s1d
          result10 = derf(s1)
        
          a1negd = -(0.5*result10d)
          a1neg = 0.5*(1-result10)
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !		Expressions for the split fluxes ..	
          gxnd = 0.D0
          gxnd(1) = rhod*(ut*a1neg-b1) + rho*(utd*a1neg+ut*a1negd-b1d)
          gxn(1) = rho*(ut*a1neg-b1)
          temp1d = pr_by_rhod + utd*ut + ut*utd
          temp1 = pr_by_rho + ut*ut
          temp2d = temp1d*a1neg + temp1*a1negd - utd*b1 - ut*b1d
          temp2 = temp1*a1neg - ut*b1
          gxnd(2) = rhod*temp2 + rho*temp2d
          gxn(2) = rho*temp2
          temp1d = (utd*un+ut*und)*a1neg + ut*un*a1negd - und*b1 - un*b1d
          temp1 = ut*un*a1neg - un*b1
          gxnd(3) = rhod*temp1 + rho*temp1d
          gxn(3) = rho*temp1
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5*((utd*temp1+ut*temp1d)*a1neg+ut*temp1*a1negd)
          temp2 = 0.5*ut*temp1*a1neg
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          gxnd(4) = rhod*(temp2-0.5*temp1*b1) + rho*(temp2d-0.5*(temp1d*b1+temp1&
        &   *b1d))
          gxn(4) = rho*(temp2-0.5*temp1*b1)
        END SUBROUTINE FLUX_GXN_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_gyp in forward (tangent) mode:
        !   variations   of useful results: gyp
        !   with respect to varying inputs: u
        !   RW status of diff variables: u:in gyp:out
        SUBROUTINE FLUX_GYP_D(gyp, gypd, u, ud, nx, ny)
          IMPLICIT NONE
          DOUBLE PRECISION :: gyp(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gypd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta, s2, b2, a2pos
          DOUBLE PRECISION :: betad, s2d, b2d, a2posd
          DOUBLE PRECISION :: temp1, temp2
          DOUBLE PRECISION :: temp1d, temp2d
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
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5*(2*u(2)*ud(2)+2*u(3)*ud(3))/u(1)+0.5*ud(1)*(u(2)**2+u&
&   (3)**2)/u(1)**2)/2.5
          pr = (u(4)-0.5/u(1)*(u(2)**2+u(3)**2))/2.5
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
          betad = (0.5*rhod*pr-0.5*rho*prd)/pr**2
          beta = 0.5*rho/pr
          IF (beta .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2posd = 0.5*result10d
          a2pos = 0.5*(1+result10)
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !		Expressions for the split fluxes ..	
          gypd = 0.D0
          gypd(1) = rhod*(un*a2pos+b2) + rho*(und*a2pos+un*a2posd+b2d)
          gyp(1) = rho*(un*a2pos+b2)
          temp1d = pr_by_rhod + und*un + un*und
          temp1 = pr_by_rho + un*un
          temp2d = temp1d*a2pos + temp1*a2posd + und*b2 + un*b2d
          temp2 = temp1*a2pos + un*b2
          gypd(3) = rhod*temp2 + rho*temp2d
          gyp(3) = rho*temp2
          temp1d = (utd*un+ut*und)*a2pos + ut*un*a2posd + utd*b2 + ut*b2d
          temp1 = ut*un*a2pos + ut*b2
          gypd(2) = rhod*temp1 + rho*temp1d
          gyp(2) = rho*temp1
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5*((und*temp1+un*temp1d)*a2pos+un*temp1*a2posd)
          temp2 = 0.5*un*temp1*a2pos
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          gypd(4) = rhod*(temp2+0.5*temp1*b2) + rho*(temp2d+0.5*(temp1d*b2+temp1&
        &   *b2d))
          gyp(4) = rho*(temp2+0.5*temp1*b2)
        END SUBROUTINE FLUX_GYP_D
        
        !        Generated by TAPENADE     (INRIA, Ecuador team)
        !  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
        !
        !  Differentiation of flux_gyn in forward (tangent) mode:
        !   variations   of useful results: gyn
        !   with respect to varying inputs: u
        !   RW status of diff variables: u:in gyn:out
        SUBROUTINE FLUX_GYN_D(gyn, gynd, u, ud, nx, ny)
          IMPLICIT NONE
          DOUBLE PRECISION :: gyn(4), u1, u2, rho, pr, u(4)
          DOUBLE PRECISION :: gynd(4), u1d, u2d, rhod, prd, ud(4)
          DOUBLE PRECISION :: tx, ty, nx, ny, ut, un
          DOUBLE PRECISION :: utd, und
          DOUBLE PRECISION :: beta, s2, b2, a2neg
          DOUBLE PRECISION :: betad, s2d, b2d, a2negd
          DOUBLE PRECISION :: temp1, temp2
          DOUBLE PRECISION :: temp1d, temp2d
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
          tx = ny
          ty = -nx
          rhod = ud(1)
          rho = u(1)
          u1d = (ud(2)*u(1)-u(2)*ud(1))/u(1)**2
          u1 = u(2)/u(1)
          u2d = (ud(3)*u(1)-u(3)*ud(1))/u(1)**2
          u2 = u(3)/u(1)
          prd = (ud(4)-0.5*(2*u(2)*ud(2)+2*u(3)*ud(3))/u(1)+0.5*ud(1)*(u(2)**2+u&
&   (3)**2)/u(1)**2)/2.5
          pr = (u(4)-0.5/u(1)*(u(2)**2+u(3)**2))/2.5
          utd = tx*u1d + ty*u2d
          ut = u1*tx + u2*ty
          und = nx*u1d + ny*u2d
          un = u1*nx + u2*ny
          betad = (0.5*rhod*pr-0.5*rho*prd)/pr**2
          beta = 0.5*rho/pr
          IF (beta .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = betad/(2.D0*DSQRT(beta))
          END IF
          result1 = DSQRT(beta)
          s2d = und*result1 + un*result1d
          s2 = un*result1
          arg1d = pi*betad
          arg1 = pi*beta
          IF (arg1 .EQ. 0.0) THEN
            result1d = 0.D0
          ELSE
            result1d = arg1d/(2.D0*DSQRT(arg1))
          END IF
          result1 = DSQRT(arg1)
          b2d = (-(0.5*(s2d*s2+s2*s2d)*DEXP(-(s2*s2))*result1)-0.5*DEXP(-(s2*s2)&
        &   )*result1d)/result1**2
          b2 = 0.5*DEXP(-(s2*s2))/result1
        !  result10d = DERF_D(s2, s2d, result10)
        
          result10d = dexp(-s2**2)*(2.d0/sqrt(pi))*s2d
          result10 = derf(s2)
        
          a2negd = -(0.5*result10d)
          a2neg = 0.5*(1-result10)
          pr_by_rhod = (prd*rho-pr*rhod)/rho**2
          pr_by_rho = pr/rho
          u_sqrd = utd*ut + ut*utd + und*un + un*und
          u_sqr = ut*ut + un*un
        !		Expressions for the split fluxes ..	
          gynd = 0.D0
          gynd(1) = rhod*(un*a2neg-b2) + rho*(und*a2neg+un*a2negd-b2d)
          gyn(1) = rho*(un*a2neg-b2)
          temp1d = pr_by_rhod + und*un + un*und
          temp1 = pr_by_rho + un*un
          temp2d = temp1d*a2neg + temp1*a2negd - und*b2 - un*b2d
          temp2 = temp1*a2neg - un*b2
          gynd(3) = rhod*temp2 + rho*temp2d
          gyn(3) = rho*temp2
          temp1d = (utd*un+ut*und)*a2neg + ut*un*a2negd - utd*b2 - ut*b2d
          temp1 = ut*un*a2neg - ut*b2
          gynd(2) = rhod*temp1 + rho*temp1d
          gyn(2) = rho*temp1
          temp1d = 7.0d0*pr_by_rhod + u_sqrd
          temp1 = 7.0d0*pr_by_rho + u_sqr
          temp2d = 0.5*((und*temp1+un*temp1d)*a2neg+un*temp1*a2negd)
          temp2 = 0.5*un*temp1*a2neg
          temp1d = 6.0d0*pr_by_rhod + u_sqrd
          temp1 = 6.0d0*pr_by_rho + u_sqr
          gynd(4) = rhod*(temp2-0.5*temp1*b2) + rho*(temp2d-0.5*(temp1d*b2+temp1&
        &   *b2d))
          gyn(4) = rho*(temp2-0.5*temp1*b2)
        END SUBROUTINE FLUX_GYN_D



end module split_fluxes_mod
