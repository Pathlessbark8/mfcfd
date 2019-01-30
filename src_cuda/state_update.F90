module state_update_mod

        use device_data_structure_mod

        contains
        
        attributes(global) subroutine state_update(x_d, nx_d, flag_d, nbhs_d, conn_d, &
                & prim_d, flux_res_d)

                implicit none

                ! device variables
                real*8 :: x_d(:,:), nx_d(:,:)
                integer :: flag_d(:,:), nbhs_d(:), conn_d(:,:)
                real*8 :: prim_d(:,:)
                real*8 :: flux_res_d(:,:)
                ! local variables
                integer :: i, k, r
                real*8 :: U(4), temp
                real*8 :: nx, ny
                real*8 :: U2_rot, U3_rot

                k = (blockIdx%x-1)* blockDim%x + threadIdx%x

                if (flag_d(1,k) == 1) then
                        
                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call primitive_to_conserved(k, nx, ny, U, prim_d)

                        U = U - flux_res_d(:,k)
                        U(3) = 0.d0

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

                if (flag_d(1,k) == 3) then

                        nx = nx_d(1,k)
                        ny = nx_d(2,k)

                        call conserved_vector_Ubar(k, U, nx, ny, prim_d) 
                        
                        U = U - flux_res_d(:,k)
                        
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

                        U = U - flux_res_d(:,k)

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

        end subroutine

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

end module state_update_mod
