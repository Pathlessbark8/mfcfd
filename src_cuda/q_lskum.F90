module q_lskum_mod
#include "limiter.h"

        use data_structure_mod
        use device_data_structure_mod
        use point_normals_mod    
        use generate_connectivity_mod
        use cudafor

contains

        subroutine q_lskum()

                implicit none

                ! Grid and block dim
                type(dim3) :: grid , tBlock
                integer :: istat
                real*8 :: residue, res_old, res_new
                real*8, device :: sum_res_sqr_d, max_res_d
                real*8 :: sum_res_sqr, max_res

                OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

                if(old_format == 0) call compute_normals()
                call generate_connectivity()

                ! Transfer from host device

                call host_to_device()

                tBlock = dim3 (blockx, blocky, blockz)
                grid = dim3(ceiling(real(max_points)/ tBlock%x), 1, 1)
                
                write(*,*)'%%%%%%%%%%%%%%%-GPU size info-%%%%%%%%%%%%%%'
                write(*,*) 'number of threads per block:',blockx*blocky*blockz
                write(*,*) 'grid dimension:',grid
                write(*,*) 'thread block dimension:',tBlock

                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-Iterations begin-%%%%%%%%%%%%%'
                write(*,*)

                do it = 1, max_iters
                        
                        max_res = 0.0d0
                        sum_res_sqr = 0.0d0
                        max_res_d = max_res
                        sum_res_sqr_d = sum_res_sqr

                        call eval_q_variables<<<grid, tBlock>>>(point_d%prim, point_d%q)

                        call eval_q_derivatives<<<grid, tBlock>>>(point_d%x, point_d%nbhs, &
                                & point_d%conn, point_d%q, point_d%dq)

                        call func_delta<<<grid, tBlock>>>(point_d%x, point_d%nbhs, &
                               & point_d%conn, point_d%prim, point_d%delta)

                        call cal_flux_residual<<<grid, tBlock>>>(point_d%x, point_d%nx, &
                                & point_d%flag, point_d%nbhs, point_d%conn, &
                                & point_d%xpos_nbhs, point_d%xneg_nbhs, point_d%ypos_nbhs, &
                                & point_d%yneg_nbhs, point_d%xpos_conn, point_d%xneg_conn,&
                                & point_d%ypos_conn, point_d%yneg_conn, point_d%prim,  &
                                & point_d%q, point_d%dq, point_d%flux_res)

                        call state_update<<<grid, tBlock>>>(point_d%x, point_d%nx, point_d%flag, &
                                & point_d%nbhs, point_d%conn, point_d%prim, point_d%flux_res, &
                                & point_d%delta, max_res_d, sum_res_sqr_d)

                        max_res = max_res_d
                        sum_res_sqr = sum_res_sqr_d

                        istat = cudaGetLastError()
                        if (istat .ne. 0) then
                                print*, cudaGetErrorString(istat)
                                stop istat
                        end if
               
                        res_new = dsqrt(sum_res_sqr)/max_points

                        if(it .le. 2) then
                                res_old = res_new
                                residue = 0.d0
                        else 
                                residue = dlog10(res_new/res_old)
                        endif

                        write(*,'(a12,i8,a15,e30.20)')'iterations:',it,'residue:',residue
                        write(301, *) it, residue

                enddo
                
                call device_to_host()

                CLOSE(UNIT=301)

        end subroutine

        attributes(global) subroutine eval_q_variables(prim_d, q_d)

                        implicit none

                        ! device variables
                        real*8 :: prim_d(:,:), q_d(:,:)
                        ! local variables
                        integer :: i
                        real*8 :: rho, u1, u2, pr, beta
                        real*8 :: two_times_beta

                        i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                        rho = prim_d(1,i)
                        u1 = prim_d(2,i)
                        u2 = prim_d(3,i)
                        pr = prim_d(4,i)

                        beta = 0.5d0*rho/pr

                        q_d(1,i) = dlog(rho) + (dlog(beta)*2.5d0) - beta*(u1*u1 + u2*u2)

                        two_times_beta = 2.0d0*beta

                        q_d(2,i) = two_times_beta*u1

                        q_d(3,i) = two_times_beta*u2

                        q_d(4,i) = -two_times_beta

                        call syncthreads()

        end subroutine

        attributes(global) subroutine eval_q_derivatives(x_d, nbhs_d, conn_d, q_d, dq_d)


                implicit none

                ! device variables
                real*8 :: q_d(:,:), dq_d(:,:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                real*8 :: x_d(:,:)
                ! local variables
                integer :: i
                integer :: k, r, nbh
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: delx, dely, dist, weights
                real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
                real*8 :: sum_delx_delq(4), sum_dely_delq(4)
                real*8 :: det, delq, temp
                real*8 :: one_by_det

                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                x_i = x_d(1,i)
                y_i = x_d(2,i)

                sum_delx_sqr = 0.d0
                sum_dely_sqr = 0.d0
                sum_delx_dely = 0.d0

                sum_delx_delq = 0.d0
                sum_dely_delq = 0.d0

                do k = 1, nbhs_d(i)

                        nbh = conn_d(i,k)

                        x_k = x_d(1,nbh)
                        y_k = x_d(2,nbh)

                        delx = x_k - x_i
                        dely = y_k - y_i

                        dist = dsqrt(delx*delx + dely*dely)
                        weights = dist**power_d

                        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                        sum_dely_sqr = sum_dely_sqr + dely*dely*weights

                        sum_delx_dely = sum_delx_dely + delx*dely*weights

                        sum_delx_delq = sum_delx_delq + weights*delx*(q_d(:,nbh) - q_d(:,i))
                        sum_dely_delq = sum_dely_delq + weights*dely*(q_d(:,nbh) - q_d(:,i))

                enddo
                call syncthreads()

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                one_by_det = 1.0d0/det

                dq_d(1,:,i) = (sum_delx_delq*sum_dely_sqr&
                       & - sum_dely_delq*sum_delx_dely)*one_by_det
                dq_d(2,:,i) = (sum_dely_delq*sum_delx_sqr&
                                       &- sum_delx_delq*sum_delx_dely)*one_by_det

                call syncthreads()

        end subroutine

        attributes(global) subroutine func_delta(x_d, nbhs_d, conn_d, prim_d, delta_d)


                implicit none
                ! device variables
                real*8 :: prim_d(:,:), delta_d(:)
                integer :: nbhs_d(:), conn_d(:,:)
                real*8 :: x_d(:,:)
                ! local variables
                integer :: i
                integer ::  k, r
                real*8 :: delta_t
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
                call syncthreads()
                delta_d(i) = min_delt

        end subroutine

        attributes(global) subroutine cal_flux_residual(x_d, nx_d, flag_d, nbhs_d, conn_d, &
                & xpos_nbhs_d, xneg_nbhs_d, ypos_nbhs_d, yneg_nbhs_d, xpos_conn_d, xneg_conn_d, &
                & ypos_conn_d, yneg_conn_d, prim_d, q_d, dq_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
                integer :: xpos_nbhs_d(:), xneg_nbhs_d(:), ypos_nbhs_d(:), yneg_nbhs_d(:)
                integer :: xpos_conn_d(:,:), xneg_conn_d(:,:), ypos_conn_d(:,:), yneg_conn_d(:,:)
                real*8 :: prim_d(:,:), q_d(:,:)
                real*8 :: flux_res_d(:,:), dq_d(:,:,:)
                ! local variables
                integer :: i
                real*8 :: Gxp(4), Gxn(4), Gyp(4), Gyn(4)

                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if (flag_d(1,i) == 1) then

                        call wall_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call wall_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call wall_dGy_neg(i, Gyn, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = Gxp + Gxn + Gyn
                        flux_res_d(:,i) = 2.0d0 * flux_res_d(:,i)
                end if
                call syncthreads()

                if (flag_d(1,i) == 3) then

                        call outer_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call outer_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call outer_dGy_pos(i, Gyp, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = Gxp + Gxn + Gyp
                end if
                call syncthreads()

                if (flag_d(1,i) == 2) then

                        call interior_dGx_pos(i, Gxp, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                        call interior_dGx_neg(i, Gxn, x_d, nx_d, nbhs_d, conn_d, xneg_nbhs_d, &
                                & xneg_conn_d, prim_d, q_d, dq_d)
                        
                        call interior_dGy_pos(i, Gyp, x_d, nx_d, nbhs_d, conn_d, ypos_nbhs_d, &
                                & ypos_conn_d, prim_d, q_d, dq_d)
                        
                        call interior_dGy_neg(i, Gyn, x_d, nx_d, nbhs_d, conn_d, yneg_nbhs_d, &
                                & yneg_conn_d, prim_d, q_d, dq_d)
                        
                        flux_res_d(:,i) = Gxp + Gxn + Gyp + Gyn
                end if
                call syncthreads()

        end subroutine
        
        attributes(global) subroutine state_update(x_d, nx_d, flag_d, nbhs_d, conn_d, &
                & prim_d, flux_res_d, delta_d, max_res, sum_res_sqr)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
                real*8 :: prim_d(:,:), delta_d(:)
                real*8 :: flux_res_d(:,:)
                real*8 :: max_res, sum_res_sqr, res_sqr
                ! local variables
                integer :: i, k, r
                real*8 :: delt, U(4), temp
                real*8 :: nx, ny
                real*8 :: U2_rot, U3_rot

                k = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if (flag_d(1,k) == 1) then
                        
                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call primitive_to_conserved(k, nx, ny, U, prim_d)

                        temp = U(1)

                        U = U - delta_d(k) * flux_res_d(:,k)
                        U(3) = 0.d0

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)
                        
                        if(res_sqr .gt. max_res) then 
                                max_res = res_sqr
                        endif

                        sum_res_sqr = sum_res_sqr + res_sqr

                        prim_d(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        prim_d(2,k) = U(2)*temp
                        prim_d(3,k) = U(3)*temp
                        prim_d(4,k) = 0.4d0*U(4) - (0.2d0*temp)*(U(2)*U(2) + U(3)*U(3))

                end if
                call syncthreads()

                if (flag_d(1,k) == 3) then

                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call conserved_vector_Ubar(k, U, nx, ny, prim_d) 
                        
                        U = U - delta_d(k) *flux_res_d(:,k)
                        
                        U2_rot = U(2)
                        
                        U3_rot = U(3)
                        
                        U(2) = U2_rot*ny + U3_rot*nx
                        
                        U(3) = U3_rot*ny - U2_rot*nx

                        prim_d(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        prim_d(2,k) = U(2)*temp
                        prim_d(3,k) = U(3)*temp
                        prim_d(4,k) = 0.4d0*U(4) - (0.2d0*temp)*(U(2)*U(2) + U(3)*U(3))


                end if
                call syncthreads()

                if (flag_d(1,k) == 2) then

                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call primitive_to_conserved(k, nx, ny, U, prim_d)

                        temp = U(1)

                        U = U - delta_d(k) * flux_res_d(:,k)

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)

                        if(res_sqr .gt. max_res) then 
                                max_res = res_sqr
                        endif

                        sum_res_sqr = sum_res_sqr + res_sqr

                        prim_d(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        prim_d(2,k) = U(2)*temp
                        prim_d(3,k) = U(3)*temp

                        prim_d(4,k) = 0.4d0*U(4) - (0.2d0*temp)*(U(2)*U(2) + U(3)*U(3))
        
                end if
                call syncthreads()

        end subroutine

        attributes(device) subroutine wall_dGx_pos(i, G, x_d, nx_d, nbhs_d, conn_d, xpos_nbhs_d, &
                                & xpos_conn_d, prim_d, q_d, dq_d)

                ! device variables
                integer :: i
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
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

                        qtilde_i = q_d(:,i) - 0.5d0*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
#ifdef VENKAT
                        call venkat_limiter(qtilde_i, phi_i, i, q_d, nbhs_d, conn_d, x_d)
                        call venkat_limiter(qtilde_k, phi_k, k, q_d, nbhs_d, conn_d, x_d)
#endif
                        qtilde_i = q_d(:,i) - 0.5d0*phi_i*(delx*dq_d(1,:,i) + dely*dq_d(2,:,i))
                        qtilde_k = q_d(:,k) - 0.5d0*phi_k*(delx*dq_d(1,:,k) + dely*dq_d(2,:,k))
                        
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

        attributes(device) subroutine primitive_to_conserved(k, nx, ny, U, prim_d)


                implicit none
                !device avriables
                real*8 :: prim_d(:,:)
                ! local variables
                real*8 :: rho
                real*8 :: U(4), nx, ny
                real*8 :: temp1, temp2
                integer :: k

                rho = prim_d(1,k)

                U(1) = rho
                temp1 = rho*prim_d(2,k)
                temp2 = rho*prim_d(3,k)
                U(4) = 2.5d0*prim_d(4,k) + 0.5d0*(temp1*temp1 + temp2*temp2)/rho

                U(2) = temp1*ny - temp2*nx
                U(3) = temp1*nx + temp2*ny


        end subroutine

        attributes(device) subroutine conserved_vector_Ubar(k, Ubar, nx, ny, prim_d)

                implicit none
                ! device variables
                real*8 :: prim_d(:,:)
                ! local variables
		real*8 :: u1_inf, u2_inf, u1_inf_rot, u2_inf_rot, e_inf
		real*8 :: u1, u2, pr, rho, u1_rot, u2_rot, e
		real*8 :: beta, S2, B2_inf, A2n_inf
		real*8 :: B2, A2p, temp1, temp2
		real*8 :: Ubar(4)
		real*8 :: nx, ny, tx, ty
                integer :: k
                real*8 :: rho_inf, pr_inf
                
                rho_inf = qinf1_d
                u1_inf = qinf2_d
                u2_inf = qinf3_d
                pr_inf = qinf4_d

                tx = ny
                ty = -nx

                u1_inf_rot = u1_inf*tx + u2_inf*ty
                u2_inf_rot = u1_inf*nx + u2_inf*ny

                temp1 = (u1_inf_rot*u1_inf_rot + u2_inf_rot*u2_inf_rot)
                e_inf = pr_inf/(rho_inf*(gamma-1.0d0)) + 0.5d0*(temp1)

                beta = (0.5d0*rho_inf)/pr_inf
                S2 = u2_inf_rot*dsqrt(beta)
                B2_inf = dexp(-S2*S2)/(2.0d0*dsqrt(pi*beta))
                A2n_inf = 0.5d0*(1.0d0-erf(S2))

                rho = prim_d(1,k)
                u1 = prim_d(2,k)
                u2 = prim_d(3,k)
                pr = prim_d(4,k)

                u1_rot = u1*tx + u2*ty
                u2_rot = u1*nx + u2*ny

                temp1 = (u1_rot*u1_rot + u2_rot*u2_rot)
                e = pr/(rho*(gamma-1.0d0)) + 0.5d0*(temp1)

                beta = (rho)/(2.0d0*pr)
                S2 = u2_rot*sqrt(beta)
                B2 = exp(-S2*S2)/(2.0d0*sqrt(pi*beta))
                A2p = 0.5d0*(1.0d0+erf(S2))

                Ubar(1) = (rho_inf*A2n_inf) + (rho*A2p)

                Ubar(2) = (rho_inf*u1_inf_rot*A2n_inf) + (rho*u1_rot*A2p)
        
                temp1 = rho_inf*(u2_inf_rot*A2n_inf - B2_inf)
                temp2 = rho*(u2_rot*A2p + B2)
                Ubar(3) = temp1 + temp2

                temp1 = (rho_inf*A2n_inf*e_inf - 0.5d0*rho_inf*u2_inf_rot*B2_inf)
                temp2 = (rho*A2p*e + 0.5d0*rho*u2_rot*B2)

                Ubar(4) = temp1 + temp2
                
        end subroutine

end module q_lskum_mod
