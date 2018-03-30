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


                do i = 1, wall_points

                        k = wall_points_index(i)

                        call wall_dGx_pos(Gxp, k) 

                        call wall_dGx_neg(Gxn, k) 

                        call wall_dGy_neg(Gyn, k) 

                        p%flux_res(:,k) = Gxp + Gxn + Gyn

                        p%flux_res(:,k) = 2.0d0*p%flux_res(:,k)



                enddo


                do i = 1, outer_points

                        k = outer_points_index(i)

                        call outer_dGx_pos(Gxp, k)

                        call outer_dGx_neg(Gxn, k) 

                        call outer_dGy_pos(Gyp, k) 

                        p%flux_res(:,k) = Gxp + Gxn + Gyp


                enddo


                do i = 1, interior_points

                        k = interior_points_index(i)

                        call interior_dGx_pos(Gxp, k) 

                        call interior_dGx_neg(Gxn, k) 

                        call interior_dGy_pos(Gyp, k) 

                        call interior_dGy_neg(Gyn, k) 

                        p%flux_res(:,k) = Gxp + Gxn + Gyp + Gyn


                enddo


        end subroutine


end module flux_residual_mod

