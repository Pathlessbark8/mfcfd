module flux_residual_mod

        use device_data_structure_mod


contains

        attributes(global) subroutine dGx_pos(x_d, nx_d, flag_d, dist_d, nbhs_d, conn_d, &
                & xpos_nbhs_d, xpos_conn_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:), dist_d(:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:)
                integer :: xpos_conn_d(:,:)
                real*8 :: q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i, j, k
                real*8 :: G_i(4), G_k(4)
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, pr, rho
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: dels_weights, deln_weights
                real*8 :: qtilde_i(4), qtilde_k(4)
                
                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if(i > mp_d) return
                
                flux_res_d(:,i) = 0.0d0

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

        !                Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d, dist_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d, dist_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        if (flag_d(i) == 0) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_quad_GxII(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_quad_GxII(G_k, nx, ny, u1, u2, rho, pr)
                        end if
                        if (flag_d(i) == 2) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_quad_GxIII(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_quad_GxIII(G_k, nx, ny, u1, u2, rho, pr)
                        end if
                        if (flag_d(i) == 1) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_Gxp(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_Gxp(G_k, nx, ny, u1, u2, rho, pr)
                        end if

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                flux_res_d(:,i) = flux_res_d(:,i) + (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det

        end subroutine
        
        attributes(global) subroutine dGx_neg(x_d, nx_d, flag_d, dist_d, nbhs_d, conn_d, &
                & xneg_nbhs_d, xneg_conn_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:), dist_d(:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                integer :: xneg_nbhs_d(:)
                integer :: xneg_conn_d(:,:)
                real*8 :: q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i, j, k
                real*8 :: G_i(4), G_k(4)
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, pr, rho
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: dels_weights, deln_weights
                real*8 :: qtilde_i(4), qtilde_k(4)
                
                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if(i > mp_d) return
                
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

        !                Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d, dist_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d, dist_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        if (flag_d(i) == 0) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_quad_GxI(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_quad_GxI(G_k, nx, ny, u1, u2, rho, pr)
                        end if
                        if (flag_d(i) == 2) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_quad_GxIV(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_quad_GxIV(G_k, nx, ny, u1, u2, rho, pr)
                        end if
                        if (flag_d(i) == 1) then
                                call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                                call flux_Gxn(G_i, nx, ny, u1, u2, rho, pr)

                                call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                                call flux_Gxn(G_k, nx, ny, u1, u2, rho, pr)
                        end if

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                flux_res_d(:,i) = flux_res_d(:,i) + (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det
                
        end subroutine
        
        attributes(global) subroutine dGy_pos(x_d, nx_d, flag_d, dist_d, nbhs_d, conn_d, &
                & ypos_nbhs_d, ypos_conn_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:), dist_d(:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                integer :: ypos_nbhs_d(:)
                integer :: ypos_conn_d(:,:)
                real*8 :: q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i, j, k, ii
                real*8 :: G_i(4), G_k(4)
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, pr, rho
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: dels_weights, deln_weights
                real*8 :: qtilde_i(4), qtilde_k(4)
                
                i = (blockIdx%x-1)* blockDim%x + threadIdx%x
                ii = threadIdx%x

                if(i > mp_d) return
                if(flag_d(i) == 0) return
                
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

        !                Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d, dist_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d, dist_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyp(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyp(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                flux_res_d(:,i) = flux_res_d(:,i) + (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det
                
        end subroutine
        
        attributes(global) subroutine dGy_neg(x_d, nx_d, flag_d, dist_d, nbhs_d, conn_d, &
                & yneg_nbhs_d, yneg_conn_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:), dist_d(:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                integer :: yneg_nbhs_d(:)
                integer :: yneg_conn_d(:,:)
                real*8 :: q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i, j, k
                real*8 :: G_i(4), G_k(4)
                real*8 :: tx, ty, nx, ny
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, pr, rho
                real*8 :: delx, dely, det, one_by_det
                real*8 :: dels, deln
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delf(4), sum_dely_delf(4)
                real*8 :: dist, weights
                real*8 :: phi_i(4), phi_k(4)
                real*8 :: dels_weights, deln_weights
                real*8 :: qtilde_i(4), qtilde_k(4)
                
                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if(i > mp_d) return
                if(flag_d(i) == 2) return
                
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

        !                Higher order accuracy using q-variables ..

                        qtilde_i = q_d(:,i) - 0.5d0*fo_flag*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*fo_flag*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d, dist_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d, dist_d)
                        
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
                        call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
                        call flux_Gyn(G_i, nx, ny, u1, u2, rho, pr)

                        call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
                        call flux_Gyn(G_k, nx, ny, u1, u2, rho, pr)

                        sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
                        sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

                end do

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.d0/det

                flux_res_d(:,i) = flux_res_d(:,i) + (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det

        end subroutine

        attributes(device) subroutine venkat_limiter(qtilde, phi, k, q_d, nbhs_d, conn_d, x_d, dist_d)


                implicit none

                ! device variables
                real*8 :: q_d(:,:), x_d(:,:), dist_d(:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: r, k
                real*8 :: qtilde(4), phi(4)
                real*8 :: q, del_neg, del_pos
                real*8 :: max_q, min_q, ds, epsi, num, den


                do r = 1, 4
                        q = q_d(r,k)
                        del_neg = qtilde(r) - q
                        if(dabs(del_neg) .le. 10e-6) then
                                phi(r)=1.d0

                        else if(dabs(del_neg) .gt. 10e-6) then
                                call max_min_q(k, r, max_q, min_q, q_d, nbhs_d, conn_d)
                                if(del_neg .gt. 0.d0) then
                                        del_pos = max_q - q
                                else if(del_neg .lt. 0.d0) then
                                        del_pos = min_q - q
                                endif

                                epsi = vl_d*dist_d(k)
                                epsi = epsi**3.0d0

                                num = (del_pos*del_pos) + (epsi*epsi)  ! Numerator .. 
                                num = num*del_neg + 2.0d0*del_neg*del_neg*del_pos

                                den = del_pos*del_pos + 2.0d0*del_neg*del_neg ! Denominator ..
                                den = den + del_neg*del_pos + epsi*epsi
                                den = den*del_neg

                                if(num/den .lt. 1.d0) then
                                        phi(r) = num/den
                                else
                                        phi(r) = 1.d0
                                endif

                        endif

                enddo

        end subroutine

        attributes(device) subroutine max_min_q(k, r, max, min, q_d, nbhs_d, conn_d)

                implicit none

                ! device variables
                real*8 :: q_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                ! local variables
                integer :: k, r, j, nbh
                real*8 :: max, min

                max = q_d(r,k)
                min = q_d(r,k)

                do j = 1, nbhs_d(k)
                        nbh = conn_d(k,j)

                        if(q_d(r,nbh) .gt. max) then
                                max = q_d(r,nbh)
                        endif

                        if(q_d(r,nbh) .lt. min) then
                                min = q_d(r,nbh)
                        endif
                enddo
        end subroutine

        attributes(device) subroutine qtilde_to_primitive(qtilde, u1, u2, rho, pr)

                implicit none

                real*8, intent(in) :: qtilde(*)
                real*8 :: u1, u2, rho, pr
                real*8 :: beta

                beta = -qtilde(4)*0.5d0

                u1 = qtilde(2)*(0.5d0/beta)
                u2 = qtilde(3)*(0.5d0/beta)

                rho = dexp((qtilde(1) + beta*(u1*u1 + u2*u2)) - (dlog(beta)/(gamma-1)))
                pr = rho*(0.5d0/beta)

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

                real*8, intent(inout) :: Gyn(*)
                real*8, intent(in) :: nx, ny, u1, u2, rho, pr
                ! local variables
                real*8 :: ut, un
                real*8 :: S2, B2, A2neg

                ut = u1*ny - u2*nx
                un = u1*nx + u2*ny

                S2 = un*dsqrt(0.5*rho/pr) 
                B2 = 0.5*dexp(-S2*S2)/dsqrt(pi*0.5*rho/pr)
                A2neg = 0.5*(1 - erf(S2))     

!		Expressions for the split fluxes ..	

                Gyn(1) = rho*(un*A2neg - B2)  

                Gyn(3) = rho*((pr/rho + un*un)*A2neg - un*B2)

                Gyn(2) = rho*(ut*un*A2neg - ut*B2)

                Gyn(4) = rho*((0.5*un*((7.0d0*pr/rho) + (ut*ut + un*un))*A2neg) - 0.5*((6.0d0*pr/rho) + (ut*ut + un*un))*B2)

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

