module flux_residual_mod


        use parameter_mod
        use data_structure_mod
        use interior_fluxes_mod
        use wall_fluxes_mod
        use outer_fluxes_mod

contains


        subroutine cal_flux_residual()


                implicit none

                integer :: i, k
		real*8 :: Gxp(4), Gxn(4), Gyp(4), Gyn(4)
		real*8 :: Gxpd(4), Gxnd(4), Gypd(4), Gynd(4)
		real*8 :: U(4)
		real*8 :: local_sos, Dp
		real*8 :: spectral_Ax, spectral_Ay
		real*8 :: Ax, Ay
                real*8 :: LSxp, LSxn, LSyp, LSyn
                do i = 1, wall_points

                        k = wall_points_index(i)

                        call wall_dGx_pos(Gxp, Gxpd, LSxp, k) 

                        call wall_dGx_neg(Gxn, Gxnd, LSxn, k) 

                        call wall_dGy_neg(Gyn, Gynd, LSyn, k)

                        local_sos = dsqrt(gamma*point%prim(4,k)/point%prim(1,k))

                        if(tscheme== 0) then
                                LSxp = 0.0d0
                                LSxn = 0.0d0
                                LSyn = 0.0d0
                                Gxpd = 0.0d0
                                Gxnd = 0.0d0
                                Gynd = 0.0d0
                        end if

                        spectral_Ax = dabs(point%prim(2,k)) + local_sos
                        spectral_Ay = dabs(point%prim(3,k)) + local_sos

                        Ax = spectral_Ax / 2.0d0
                        Ay = spectral_Ay / 2.0d0

                        Dp = 1/point%delta(k) - (LSxp - LSxn) * Ax - (-LSyn) * Ay

                        point%flux_res(:,k) = 1/Dp*(Gxp + Gxn + Gyn + Gxpd + Gxnd + Gynd)
                        
                        point%flux_res(:,k) = 2.0d0*point%flux_res(:,k)

                enddo

                do i = 1, outer_points

                        k = outer_points_index(i)

                        call outer_dGx_pos(Gxp, Gxpd, LSxp,  k)

                        call outer_dGx_neg(Gxn, Gxnd, LSxn, k) 

                        call outer_dGy_pos(Gyp, Gypd, LSyp, k) 

                        local_sos = dsqrt(gamma*point%prim(4,k)/point%prim(1,k))
                        
                        if(tscheme== 0) then
                                LSxp = 0.0d0
                                LSxn = 0.0d0
                                LSyp = 0.0d0
                                Gxpd = 0.0d0
                                Gxnd = 0.0d0
                                Gypd = 0.0d0
                        end if
                        
                        spectral_Ax = dabs(point%prim(2,k)) + local_sos
                        spectral_Ay = dabs(point%prim(3,k)) + local_sos

                        Ax = spectral_Ax / 2.0d0
                        Ay = spectral_Ay / 2.0d0
                        
                        Dp = 1/point%delta(k) - (LSxp - LSxn) * Ax - (LSyp) * Ay

                        point%flux_res(:,k) = 1/Dp*(Gxp + Gxn + Gyp + Gxpd + Gxnd + Gypd)
                enddo


                do i = 1, interior_points

                        k = interior_points_index(i)
                        
                        LSxp = 0.0d0
                        LSxn = 0.0d0
                        LSyp = 0.0d0
                        LSyn = 0.0d0
                        Gxpd = 0.0d0
                        Gxnd = 0.0d0
                        Gypd = 0.0d0
                        Gynd = 0.0d0

                        call interior_dGx_pos(Gxp, Gxpd, LSxp, k) 

                        call interior_dGx_neg(Gxn, Gxnd, LSxn, k) 

                        call interior_dGy_pos(Gyp, Gypd, LSyp, k) 

                        call interior_dGy_neg(Gyn, Gynd, LSyn, k) 

                        local_sos = dsqrt(gamma*point%prim(4,k)/point%prim(1,k))
                        
                        if(tscheme== 0) then
                                LSxp = 0.0d0
                                LSxn = 0.0d0
                                LSyp = 0.0d0
                                LSyn = 0.0d0
                                Gxpd = 0.0d0
                                Gxnd = 0.0d0
                                Gypd = 0.0d0
                                Gynd = 0.0d0
                        end if
                        
                        spectral_Ax = dabs(point%prim(2,k)) + local_sos
                        spectral_Ay = dabs(point%prim(3,k)) + local_sos

                        Ax = spectral_Ax / 2.0d0
                        Ay = spectral_Ay / 2.0d0

                        Dp = 1/point%delta(k) - (LSxp - LSxn) * Ax - (LSyp - LSyn) * Ay

                        point%flux_res(:,k) = 1/Dp*(Gxp + Gxn + Gyp + Gyn + Gxpd + Gxnd + Gypd + Gynd)
                
                enddo

                ! Save previous solution
                do i=1,local_points
                        point%U_old(1,i) = point%prim(1,i)
                        point%U_old(2,i) = point%prim(1,i)*point%prim(2,i)
                        point%U_old(3,i) = point%prim(1,i)*point%prim(3,i)
                        point%U_old(4,i) = 2.5d0*point%prim(4,i) + 0.5d0*point%prim(1,i)*&
                                &(point%prim(2,i)*point%prim(2,i) +&
                                &point%prim(3,i)*point%prim(3,i))
                end do

        end subroutine

end module flux_residual_mod

