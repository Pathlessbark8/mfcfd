module state_update_mod



        use data_structure_mod
        use flux_residual_mod


        contains

        
        subroutine state_update()

                implicit none

                integer :: i, k, r
		real*8 :: delt, U(4), temp
		real*8 :: res_sqr
		real*8 :: nx, ny
		real*8 :: U2_rot, U3_rot

                max_res = 0.0d0	
                sum_res_sqr = 0.0d0


                do i = 1, wall_points

                        k = wall_points_index(i)

                        nx = point%nx(k)
                        ny = point%ny(k)

                        call primitive_to_conserved(k, nx, ny, U)



                        temp = U(1)

                        U = U - point%delta(k)*point%flux_res(:,k)
                        U(3) = 0.d0

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)
                        
                        if(res_sqr .gt. max_res) then 
                                max_res = res_sqr
                                max_res_point = k
                        endif

                        sum_res_sqr = sum_res_sqr + res_sqr

                        point%prim(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        point%prim(2,k) = U(2)*temp
                        point%prim(3,k) = U(3)*temp
                        point%prim(4,k) = 0.4*U(4) - (0.2*temp)*(U(2)*U(2) + U(3)*U(3))

                enddo

                do i = 1, outer_points

                        k = outer_points_index(i)

                        nx = point%nx(k)
                        ny = point%ny(k)

                        call conserved_vector_Ubar(k, U, nx, ny) 
                        
                        temp = U(1)

                        U = U - point%delta(k)*point%flux_res(:,k)
                        
                        U2_rot = U(2)
                        
                        U3_rot = U(3)
                        
                        U(2) = U2_rot*ny + U3_rot*nx
                        
                        U(3) = U3_rot*ny - U2_rot*nx


!	Note: The outer points are not contributing to the 
!	residue ..	

!							res_sqr = (U(1) - temp)*(U(1) - temp)					
!							if(res_sqr .gt. max_res) then 
!								max_res = res_sqr
!								max_res_point = k
!							endif
!							
!							sum_res_sqr = sum_res_sqr + res_sqr

                        point%prim(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        point%prim(2,k) = U(2)*temp
                        point%prim(3,k) = U(3)*temp
                        point%prim(4,k) = 0.4*U(4) - (0.2*temp)*(U(2)*U(2) + U(3)*U(3))

                enddo

                do i = 1, interior_points
                
                        k = interior_points_index(i)

                        nx = point%nx(k)
                        ny = point%ny(k)

                        call primitive_to_conserved(k, nx, ny, U)


                        temp = U(1)
                        U = U - point%delta(k)*point%flux_res(:,k)

                        U2_rot = U(2)
                        U3_rot = U(3)
                        U(2) = U2_rot*ny + U3_rot*nx
                        U(3) = U3_rot*ny - U2_rot*nx

                        res_sqr = (U(1) - temp)*(U(1) - temp)

                        if(res_sqr .gt. max_res) then 
                                max_res = res_sqr
                                max_res_point = k
                        endif

                        sum_res_sqr = sum_res_sqr + res_sqr

                        point%prim(1,k) = U(1)
                        temp = 1.0d0/U(1)
                        point%prim(2,k) = U(2)*temp
                        point%prim(3,k) = U(3)*temp

                        point%prim(4,k) = 0.4*U(4) - (0.2*temp)*(U(2)*U(2) + U(3)*U(3))
        
                enddo




        end subroutine


        subroutine primitive_to_conserved(k, nx, ny, U) 


                implicit none
                
                
		real*8 :: rho
		real*8 :: U(4), nx, ny
		real*8 :: temp1, temp2
                integer :: k

                rho = point%prim(1,k)

                U(1) = rho 
                temp1 = rho*point%prim(2,k)
                temp2 = rho*point%prim(3,k)
                U(4) = 2.5*point%prim(4,k) + 0.5*(temp1*temp1 + temp2*temp2)/rho

                U(2) = temp1*ny - temp2*nx
                U(3) = temp1*nx + temp2*ny


        end subroutine



        subroutine conserved_to_primitive(U, k) 

                implicit none
        
		real*8 :: temp, U(4)
                integer :: k

                point%prim(1,k) = U(1)

                temp = 1.0d0/U(1)

                point%prim(2,k) = U(2)*temp
                point%prim(3,k) = U(3)*temp

                temp = U(4) - (0.5*temp)*(U(2)*U(2) + U(3)*U(3))

                point%prim(4,k) = 0.4*temp


        end subroutine


!	This subroutine computes the delta_t (local time step) at a given point ..


        subroutine func_delta()


                implicit none

                integer :: i, k, r
		real*8 :: delta_t
		real*8 :: min_dist
		real*8 :: x_i, y_i, x_k, y_k
		real*8 :: u1, u2, rho, pr, mod_u
		real*8 :: dist
		real*8 :: min_delt 

                do i = 1,local_points
                        min_delt = 1.0d0


                        do r = 1, point%nbhs(i)
                                k = point%conn(i,r)

                                rho = point%prim(1,k)
                                u1 = point%prim(2,k)
                                u2 = point%prim(3,k)
                                pr = point%prim(4,k)

                                x_i = point%x(i)
                                y_i = point%y(i)

                                x_k = point%x(k)
                                y_k = point%y(k)

                                dist = (x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i)
                                dist = dsqrt(dist)

                                mod_u = sqrt(u1*u1 + u2*u2)

                                delta_t = dist/(mod_u + 3*dsqrt(pr/rho))

                                delta_t = cfl*delta_t

                                if(min_delt > delta_t) then 
                                        min_delt = delta_t
                                endif

                        enddo
                        point%delta(i) = min_delt
                end do


        end subroutine


        subroutine conserved_vector_Ubar(k, Ubar, nx, ny)

                implicit none


		real*8 :: u1_inf, u2_inf, u1_inf_rot, u2_inf_rot, e_inf
		real*8 :: u1, u2, pr, rho, u1_rot, u2_rot, e
		real*8 :: beta, S2, B2_inf, A2n_inf
		real*8 :: B2, A2p, temp1, temp2
		real*8 :: Ubar(4)
		real*8 :: nx, ny, tx, ty
                integer :: k

                u1_inf = Mach*dcos(theta)
                u2_inf = Mach*dsin(theta)

                tx = ny
                ty = -nx

                u1_inf_rot = u1_inf*tx + u2_inf*ty
                u2_inf_rot = u1_inf*nx + u2_inf*ny

                temp1 = (u1_inf_rot*u1_inf_rot + u2_inf_rot*u2_inf_rot)
                e_inf = pr_inf/(rho_inf*(gamma-1)) + 0.5*(temp1)

                beta = (0.5*rho_inf)/pr_inf
                S2 = u2_inf_rot*dsqrt(beta)
                B2_inf = exp(-S2*S2)/(2*dsqrt(pi*beta))
                A2n_inf = 0.5*(1-derf(S2))

                rho = point%prim(1,k)
                u1 = point%prim(2,k)
                u2 = point%prim(3,k)
                pr = point%prim(4,k)

                u1_rot = u1*tx + u2*ty
                u2_rot = u1*nx + u2*ny

                temp1 = (u1_rot*u1_rot + u2_rot*u2_rot)
                e = pr/(rho*(gamma-1)) + 0.5*(temp1)

                beta = (rho)/(2*pr)
                S2 = u2_rot*sqrt(beta)
                B2 = exp(-S2*S2)/(2*sqrt(pi*beta))
                A2p = 0.5*(1+derf(S2))

                Ubar(1) = (rho_inf*A2n_inf) + (rho*A2p)

                Ubar(2) = (rho_inf*u1_inf_rot*A2n_inf) + (rho*u1_rot*A2p)
        
                temp1 = rho_inf*(u2_inf_rot*A2n_inf - B2_inf)
                temp2 = rho*(u2_rot*A2p + B2)
                Ubar(3) = temp1 + temp2

                temp1 = (rho_inf*A2n_inf*e_inf - 0.5*rho_inf*u2_inf_rot*B2_inf)
                temp2 = (rho*A2p*e + 0.5*rho*u2_rot*B2)

                Ubar(4) = temp1 + temp2
                

        end subroutine


end module state_update_mod
