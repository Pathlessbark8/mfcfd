module state_update_mod

        use device_data_structure_mod

        contains
        
        attributes(global) subroutine state_update(x_d, nx_d, flag_d, nbhs_d, conn_d, &
                & prim_d, primold_d, delta_d, flux_res_d, sum_res_sqr, rk)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:), nbhs_d(:), conn_d(:,:)
                real*8 :: prim_d(:,:), primold_d(:,:), delta_d(:)
                real*8 :: flux_res_d(:,:), sum_res_sqr(:)
                integer, value :: rk
                ! local variables
                integer :: i, k, r
                real*8 :: U(4), temp, Uold(4)
                real*8 :: nx, ny
                real*8 :: U2_rot, U3_rot, res_sqr
                integer :: istat
                real*8, parameter :: obt = 1.0d0/3.0d0
                real*8, parameter :: tbt = 2.0d0/3.0d0

                k = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if(k > mp_d) return

                sum_res_sqr(k) = 0.0d0
                
                if (flag_d(k) == 0) then
                        
                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call primitive_to_conserved(nx, ny, U, prim_d(:,k))
                        call primitive_to_conserved(nx, ny, Uold, primold_d(:,k))

                        temp = U(1)

                        if(rk == 1 .or. rk == 2 .or. rk == 4) then
                                U = U - 0.5d0 * eu_d * delta_d(k) * 2.0d0 * flux_res_d(:,k)
                        elseif(rk == 3) then
                                U = tbt * Uold + obt * (U - 0.5d0 * delta_d(k) * 2.0d0 * flux_res_d(:,k))
                        end if

                        U(3) = 0.d0

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)

                        sum_res_sqr(k) = res_sqr
                        
                        prim_d(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        prim_d(2,k) = U(2)*temp
                        prim_d(3,k) = U(3)*temp
                        prim_d(4,k) = 0.4d0*U(4) - (0.2d0*temp)*(U(2)*U(2) + U(3)*U(3))

                end if

                if (flag_d(k) == 2) then

                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call conserved_vector_Ubar(U, nx, ny, prim_d(:,k)) 
                        call conserved_vector_Ubar(Uold, nx, ny, primold_d(:,k)) 
                        
                        if(rk == 1 .or. rk == 2 .or. rk == 4) then
                                U = U - 0.5d0 * eu_d * delta_d(k) * flux_res_d(:,k)
                        elseif(rk == 3) then
                                U = tbt * Uold + obt * (U - 0.5d0 * delta_d(k) * flux_res_d(:,k))
                        end if
                        
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

                if (flag_d(k) == 1) then

                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call primitive_to_conserved(nx, ny, U, prim_d(:,k))
                        call primitive_to_conserved(nx, ny, Uold, primold_d(:,k))

                        temp = U(1)

                        if(rk == 1 .or. rk == 2 .or. rk == 4) then
                                U = U - 0.5d0 * eu_d * delta_d(k) * flux_res_d(:,k)
                        elseif(rk == 3) then
                                U = tbt * Uold + obt * (U - 0.5d0 * delta_d(k) * flux_res_d(:,k))
                        end if

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)

                        sum_res_sqr(k) = res_sqr

                        prim_d(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        prim_d(2,k) = U(2)*temp
                        prim_d(3,k) = U(3)*temp

                        prim_d(4,k) = 0.4d0*U(4) - (0.2d0*temp)*(U(2)*U(2) + U(3)*U(3))
        
                end if

        end subroutine

        attributes(device) subroutine primitive_to_conserved( nx, ny, U, prim_d)


                implicit none
                ! local variables
                real*8 :: rho, prim_d(4)
                real*8 :: U(4), nx, ny
                real*8 :: temp1, temp2

                rho = prim_d(1)

                U(1) = rho
                temp1 = rho*prim_d(2)
                temp2 = rho*prim_d(3)
                U(4) = 2.5d0*prim_d(4) + 0.5d0*(temp1*temp1 + temp2*temp2)/rho

                U(2) = temp1*ny - temp2*nx
                U(3) = temp1*nx + temp2*ny


        end subroutine

        attributes(device) subroutine conserved_vector_Ubar(Ubar, nx, ny, prim_d)

                implicit none
                ! local variables
		real*8 :: u1_inf, u2_inf, u1_inf_rot, u2_inf_rot, e_inf
		real*8 :: u1, u2, pr, rho, u1_rot, u2_rot, e, prim_d(4)
		real*8 :: beta, S2, B2_inf, A2n_inf
		real*8 :: B2, A2p, temp1, temp2
		real*8 :: Ubar(4)
		real*8 :: nx, ny, tx, ty
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

                rho = prim_d(1)
                u1 = prim_d(2)
                u2 = prim_d(3)
                pr = prim_d(4)

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

        attributes(global) subroutine eval_timestep(x_d, nbhs_d, conn_d, delta_d, prim_d, primold_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:)
                integer :: nbhs_d(:), conn_d(:,:)
                real*8 :: prim_d(:,:), delta_d(:)
                real*8 :: primold_d(:,:)
                ! local variables
                integer :: i, r, k
                real*8 :: delta_t
                real*8 :: min_dist
                real*8 :: x_i, y_i, x_k, y_k
                real*8 :: u1, u2, rho, pr, mod_u
                real*8 :: dist
                real*8 :: min_delt

                i = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if(i > mp_d) return
                
                primold_d(:,i) = prim_d(:,i)
                
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
                delta_d(i) = min_delt

        end subroutine






end module state_update_mod
