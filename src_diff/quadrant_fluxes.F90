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
                double precision derf


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
                double precision derf


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
                double precision derf


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
                double precision derf
           
          
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


end module quadrant_fluxes_mod
