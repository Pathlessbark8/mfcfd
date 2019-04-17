module flux_residual_mod
#include "limiter.h"

        use device_data_structure_mod


contains

        attributes(global) subroutine cal_flux_residual(x_d, nx_d, flag_d, nbhs_d, conn_d, &
                & xpos_nbhs_d, xneg_nbhs_d, ypos_nbhs_d, yneg_nbhs_d, xpos_conn_d, xneg_conn_d, &
                & ypos_conn_d, yneg_conn_d, prim_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:), xneg_nbhs_d(:), ypos_nbhs_d(:), yneg_nbhs_d(:)
                integer :: xpos_conn_d(:,:), xneg_conn_d(:,:), ypos_conn_d(:,:), yneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i, r, k
                real*8 :: Gxp(4), Gxn(4), Gyp(4), Gyn(4)
                ! delta t variables
                real*8 :: delta_t, delta_d
                real*8 :: min_dist
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, rho, pr, mod_u
                real*8 :: dist
                real*8 :: min_delt

                i = (blockIdx%x-1)* blockDim%x + threadIdx%x
                
                min_delt = 1.0d0

                do r = 1, nbhs_d(i)
                        k = conn_d(i,r)

                        rho = prim_d(1,k)
                        u1 = prim_d(2,k)
                        u2 = prim_d(3,k)
                        pr = prim_d(4,k)

                        x_i = x_d(1,i)
                        y_i = x_d(2,i)

                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        dist = (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i)
                        dist = dsqrt(dist)

                        mod_u = dsqrt(u1*u1 + u2*u2)

                        delta_t = dist/(mod_u + 3.0d0*dsqrt(pr/rho))

                        delta_t = cfl_d*delta_t

                        if(min_delt > delta_t) then
                                min_delt = delta_t
                        endif

                enddo
                delta_d = min_delt
                call syncthreads()

                if (flag_d(i) == 0) then

                        call wall_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call wall_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call wall_dGy_neg(i, Gyn, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = Gxp + Gxn + Gyn
                        flux_res_d(:,i) = 2.0d0 * delta_d * flux_res_d(:,i)
                end if
                call syncthreads()

                if (flag_d(i) == 2) then

                        call outer_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call outer_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call outer_dGy_pos(i, Gyp, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = delta_d * (Gxp + Gxn + Gyp)
                end if
                call syncthreads()

                if (flag_d(i) == 1) then

                        call interior_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call interior_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call interior_dGy_pos(i, Gyp, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)
                        
                        call interior_dGy_neg(i, Gyn, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = delta_d * (Gxp + Gxn + Gyp + Gyn)
                end if
                call syncthreads()

        end subroutine
        
        attributes(device) subroutine wall_dGx_pos(i, G, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:)
                integer :: xpos_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xpos_nbhs_d(i)

                        k = xpos_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif

                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_quad_GxII(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_quad_GxII(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

        attributes(device) subroutine wall_dGx_neg(i, G, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xneg_nbhs_d(:)
                integer :: xneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xneg_nbhs_d(i)

                        k = xneg_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_quad_GxI(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_quad_GxI(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

        attributes(device) subroutine wall_dGy_neg(i, G, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: yneg_nbhs_d(:)
                integer :: yneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, yneg_nbhs_d(i)

                        k = yneg_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyn(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyn(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

        attributes(device) subroutine outer_dGx_pos(i, G, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:)
                integer :: xpos_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xpos_nbhs_d(i)

                        k = xpos_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_quad_GxIII(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_quad_GxIII(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

        attributes(device) subroutine outer_dGx_neg(i, G, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xneg_nbhs_d(:)
                integer :: xneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xneg_nbhs_d(i)

                        k = xneg_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_quad_GxIV(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_quad_GxIV(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine
        
        attributes(device) subroutine outer_dGy_pos(i, G, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: ypos_nbhs_d(:)
                integer :: ypos_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, ypos_nbhs_d(i)

                        k = ypos_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyp(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyp(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine
        
        attributes(device) subroutine interior_dGx_pos(i, G, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:)
                integer :: xpos_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xpos_nbhs_d(i)

                        k = xpos_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gxp(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gxp(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

        attributes(device) subroutine interior_dGx_neg(i, G, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: xneg_nbhs_d(:)
                integer :: xneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, xneg_nbhs_d(i)

                        k = xneg_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gxn(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gxn(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine
        
        attributes(device) subroutine interior_dGy_pos(i, G, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: ypos_nbhs_d(:)
                integer :: ypos_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, ypos_nbhs_d(i)

                        k = ypos_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyp(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyp(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine
        
        attributes(device) subroutine interior_dGy_neg(i, G, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                integer :: yneg_nbhs_d(:)
                integer :: yneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: j, k, r
                real*8 :: rho, u1, u2, pr
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: G_i(4), G_k(4), G(4)
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: temp, qtilde_i(4), qtilde_k(4)
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: maxi(4), mini(4)
                real*8 :: dels_weights, deln_weights

                sum_delx_sqr = 0.0d0
                sum_dely_sqr = 0.0d0
                sum_delx_dely = 0.0d0

                sum_delx_delf = 0.0d0
                sum_dely_delf = 0.0d0

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                nx = nx_d(1,i)
                ny = nx_d(2,i)

                tx = ny
                ty = -nx

                do j = 1, yneg_nbhs_d(i)

                        k = yneg_conn_d(i,j)


                        x_k = x_d(1,k)
                        y_k = x_d(2,k)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny

                        dist = dsqrt(dels*dels + deln*deln)
                        weights = dist**power_d

                        dels_weights = dels*weights
                        deln_weights = deln*weights

                        sum_delx_sqr = sum_delx_sqr + dels*dels_weights
                        sum_dely_sqr = sum_dely_sqr + deln*deln_weights

                        sum_delx_dely = sum_delx_dely + dels*deln_weights

!                       Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#endif
#ifdef MINMAX
                        call max_q_value(i, maxi, q_d, nbhs_d, conn_d)
                        call min_q_value(i, mini, q_d, nbhs_d, conn_d)

                        do r = 1, 4
                                if( qtilde_i(r) .gt. maxi(r) ) then
                                        qtilde_i(r) = maxi(r)
                                endif

                                if( qtilde_i(r) .lt. mini(r) ) then
                                        qtilde_i(r) = mini(r)
                                endif

                                if( qtilde_k(r) .gt. maxi(r) ) then
                                        qtilde_k(r) = maxi(r)
                                endif

                                if( qtilde_k(r) .lt. mini(r) ) then
                                        qtilde_k(r) = mini(r)
                                endif
                        enddo
                        
#endif
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyn(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyn(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det
                call syncthreads()

        end subroutine

#ifdef VENKAT

        attributes(device) subroutine venkat_limiter(qtilde, phi, k, q_d, nbhs_d, conn_d, x_d)


                implicit none

                ! device variables
                real*8 :: q_d(:,:), x_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: r, k
                real*8 :: qtilde(4), phi(4)
                real*8 :: q, del_neg, del_pos
                real*8 :: max_q, min_q, ds, epsi, num, den, temp

                do r = 1, 4

                        q = q_d(r,k)
                        del_neg = qtilde(r) - q
                        if(dabs(del_neg) .le. 10e-6) then
                                phi(r)=1.d0

                        else if(dabs(del_neg) .gt. 10e-6) then
                                if(del_neg .gt. 0.d0) then
                                        call maximum(k, r, max_q, q_d, nbhs_d, conn_d)
                                        del_pos = max_q - q
                                else if(del_neg .lt. 0.d0) then
                                        call minimum(k, r, min_q, q_d, nbhs_d, conn_d)
                                        del_pos = min_q - q
                                endif

                                call smallest_dist(k, ds, x_d, nbhs_d, conn_d)

                                epsi = vl_d*ds
                                epsi = epsi**3.0d0

                                num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
                                num = num*del_neg + 2.0d0*del_neg*del_neg*del_pos

                                den = del_pos*del_pos + 2.0d0*del_neg*del_neg ! Denominator ..
                                den = den + del_neg*del_pos + epsi*epsi
                                den = den*del_neg

                                temp = num/den

                                if(temp .lt. 1.d0) then
                                        phi(r) = temp
                                else
                                        phi(r) = 1.d0
                                endif

                        endif

                enddo
                call syncthreads()

        end subroutine

        attributes(device) subroutine maximum(k, r, max, q_d, nbhs_d, conn_d)

                implicit none

                ! device variables
                real*8 :: q_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: k, r, j, nbh
                real*8 :: max

                max = q_d(r,k)

                do j = 1, nbhs_d(k)
                        nbh = conn_d(k,j)

                        if(q_d(r,nbh) .gt. max) then
                                max = q_d(r,nbh)
                        endif
                enddo
                call syncthreads()
        end subroutine

        attributes(device) subroutine minimum(k, r, min, q_d, nbhs_d, conn_d)

                implicit none

                ! device variables
                real*8 :: q_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: k, r, j, nbh
                real*8 :: min

                min = q_d(r,k)

                do j = 1, nbhs_d(k)
                        nbh = conn_d(k,j)

                        if(q_d(r,nbh) .lt. min) then
                                min = q_d(r,nbh)
                        endif
                enddo
                call syncthreads()
        end subroutine

        attributes(device) subroutine smallest_dist(k, min_dist, x_d, nbhs_d, conn_d)

                implicit none
                ! device variables
                real*8 :: x_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: k, j, nbh
                real*8 :: dx, dy, ds, min_dist

                min_dist = 10000.d0

                do j = 1, nbhs_d(k)
                        nbh = conn_d(k,j)
                        dx = x_d(1,nbh) - x_d(1,k)
                        dy = x_d(2,nbh) - x_d(2,k)

                        ds = dsqrt(dx*dx + dy*dy)

                        if(ds .lt. min_dist) then
                                min_dist = ds
                         endif

                enddo
                call syncthreads()

        end subroutine

#endif

#ifdef MINMAX

        attributes(device) subroutine max_q_value(i, maxi, q_d, nbhs_d, conn_d)

                implicit none
                ! device variables
                real*8 :: q_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer i, j, k, r
                real*8 :: maxi(4)

                maxi = q_d(:,i)

                do j = 1, nbhs_d(i)

                        k = conn_d(i,j)

                        do r = 1, 4
                                if( maxi(r) < q_d(r,k) ) then
                                        maxi(r) = q_d(r,k)
                                endif
                        enddo
                enddo

        end subroutine

        attributes(device) subroutine min_q_value(i, mini, q_d, nbhs_d, conn_d)

                implicit none
                ! device variables
                real*8 :: q_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer i, j, k, r
                real*8 :: mini(4)

                mini = q_d(:,i)

                do j = 1, nbhs_d(i)

                        k = conn_d(i,j)

                        do r = 1, 4
                                if( mini(r) > q_d(r,k) ) then
                                        mini(r) = q_d(r,k)
                                endif
                        enddo
                enddo

        end subroutine

#endif

        attributes(device) subroutine qtilde_to_primitive(qtilde, u1, u2, rho, pr)

                implicit none

                real*8 :: qtilde(4), u1, u2, rho, pr
                real*8 :: beta, temp, temp1, temp2
                real*8 :: q1, q2, q3, q4

                q1 = qtilde(1)
                q2 = qtilde(2)
                q3 = qtilde(3)
                q4 = qtilde(4)

                beta = -q4*0.5d0

                temp = 0.5d0/beta

                u1 = q2*temp
                u2 = q3*temp

                temp1 = q1 + beta*(u1*u1 + u2*u2)
                temp2 = temp1 - (dlog(beta)/(gamma-1))

                rho = dexp(temp2)
                pr = rho*temp

        end subroutine

        attributes(device) subroutine flux_Gxp(Gxp, nx, ny, u1, u2, rho, pr)

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
                A1pos = 0.5*(1 + erf(S1))     

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


        attributes(device) subroutine flux_Gxn(Gxn, nx, ny, u1, u2, rho, pr)


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
                A1neg = 0.5*(1 - erf(S1))     

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


        attributes(device) subroutine flux_Gyp(Gyp, nx, ny, u1, u2, rho, pr)


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
                A2pos = 0.5*(1 + erf(S2))     

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


        attributes(device) subroutine flux_Gyn(Gyn, nx, ny, u1, u2, rho, pr)


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
                A2neg = 0.5*(1 - erf(S2))     

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

        end subroutine

        attributes(device) subroutine flux_quad_GxI(G, nx, ny, u1, u2, rho, pr)


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
                A1neg = 0.5d0*(1.0d0 - erf(S1))     
                A2neg = 0.5d0*(1.0d0 - erf(S2))     

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



        attributes(device) subroutine flux_quad_GxII(G, nx, ny, u1, u2, rho, pr)


                implicit none

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
                A1pos = 0.5d0*(1.d0 + erf(S1))     
                A2neg = 0.5d0*(1.d0 - erf(S2))     

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



        attributes(device) subroutine flux_quad_GxIII(G, nx, ny, u1, u2, rho, pr)


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
                A1pos = 0.5d0*(1.0d0 + erf(S1))     
                A2pos = 0.5d0*(1.0d0 + erf(S2))     

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



        attributes(device) subroutine flux_quad_GxIV(G, nx, ny, u1, u2, rho, pr)


                implicit none

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
                A1neg = 0.5d0*(1.0d0 - erf(S1))     
                A2pos = 0.5d0*(1.0d0 + erf(S2))     
          
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

end module flux_residual_mod

