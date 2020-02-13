module interior_fluxes_mod

!	First written on 14.10.2016.
!	updated on Sep 23, 2017
!	

    use data_structure_mod
    use split_fluxes_mod    
    use q_variables_mod
    use limiters_mod


contains


!	This subroutine evaluates the interior flux derivative dGx_pos


    subroutine interior_dGx_pos(G, i)


        implicit none

        integer :: i, j, k, r
		real*8 :: rho, u1, u2, pr
		real*8 :: tx, ty, nx, ny
		real*8 :: x_i, y_i, x_k, y_k
		real*8 :: G_i(4), G_k(4), G(4)
		real*8 :: delx, dely, det, one_by_det
		real*8 :: dels, deln
!		
		real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
		real*8 :: sum_delx_delf(4), sum_dely_delf(4)
		real*8 :: dist, weights
		real*8 :: temp, qtilde_i(4), qtilde_k(4)
		real*8 :: phi_i(4), phi_k(4)
		real*8 :: dels_weights, deln_weights
		real*8 :: maxi(4), mini(4)


        sum_delx_sqr = 0.0d0
        sum_dely_sqr = 0.0d0
        sum_delx_dely = 0.0d0

        sum_delx_delf = 0.0d0
        sum_dely_delf = 0.0d0

        x_i = point%x(i)
        y_i = point%y(i)

        nx = point%nx(i)
        ny = point%ny(i)

        tx = ny
        ty = -nx

        do j = 1, point%xpos_nbhs(i)

            k = point%xpos_conn(i,j)

            x_k = point%x(k)
            y_k = point%y(k)

            delx = x_k - x_i
            dely = y_k - y_i

            dels = delx*tx + dely*ty
            deln = delx*nx + dely*ny

            dist = dsqrt(dels*dels + deln*deln)
            weights = dist**power

            dels_weights = dels*weights
            deln_weights = deln*weights

            sum_delx_sqr = sum_delx_sqr + dels*dels_weights
            sum_dely_sqr = sum_dely_sqr + deln*deln_weights

            sum_delx_dely = sum_delx_dely + dels*deln_weights


            qtilde_i = point%q(:,i) - 0.5d0*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))

            call venkat_limiter(qtilde_i, phi_i, i)
            call venkat_limiter(qtilde_k, phi_k, k)

            qtilde_i = point%q(:,i) - 0.5d0*phi_i*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*phi_k*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))


            call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
            call flux_Gxp(G_i, nx, ny, u1, u2, rho, pr)

            call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
            call flux_Gxp(G_k, nx, ny, u1, u2, rho, pr)
    
            sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
            sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

        enddo


        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
        one_by_det = 1.d0/det

        G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det


    end subroutine



!	This subroutine evaluates the interior flux derivative dGx_neg


    subroutine interior_dGx_neg(G, i)


        implicit none

        integer :: i, j, k, r
		real*8 :: rho, u1, u2, pr
		real*8 :: x_i, y_i, x_k, y_k
		real*8 :: tx, ty, nx, ny
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

        x_i = point%x(i)
        y_i = point%y(i)


        nx = point%nx(i)
        ny = point%ny(i)

        tx = ny
        ty = -nx

        do j = 1, point%xneg_nbhs(i)

            k = point%xneg_conn(i,j)

            x_k = point%x(k)
            y_k = point%y(k)

            delx = x_k - x_i
            dely = y_k - y_i

            dels = delx*tx + dely*ty
            deln = delx*nx + dely*ny

            dist = dsqrt(dels*dels + deln*deln)
            weights = dist**power

            dels_weights = dels*weights
            deln_weights = deln*weights

            sum_delx_sqr = sum_delx_sqr + dels*dels_weights
            sum_dely_sqr = sum_dely_sqr + deln*deln_weights

            sum_delx_dely = sum_delx_dely + dels*deln_weights


            qtilde_i = point%q(:,i) - 0.5d0*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))

            call venkat_limiter(qtilde_i, phi_i, i)
            call venkat_limiter(qtilde_k, phi_k, k)

            qtilde_i = point%q(:,i) - 0.5d0*phi_i*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*phi_k*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))


            call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
            call flux_Gxn(G_i, nx, ny, u1, u2, rho, pr)

            call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
            call flux_Gxn(G_k, nx, ny, u1, u2, rho, pr)

            sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
            sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

        enddo

        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
        one_by_det = 1.d0/det

        G = (sum_delx_delf*sum_dely_sqr - sum_dely_delf*sum_delx_dely)*one_by_det


        
    end subroutine



!	This subroutine evaluates the interior flux derivative dGx_neg
!
!
    subroutine interior_dGy_pos(G, i)


        implicit none

        integer :: i, j, k, r
		real*8 :: rho, u1, u2, pr
		real*8 :: x_i, y_i, x_k, y_k
		real*8 :: tx, ty, nx, ny

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

        x_i = point%x(i)
        y_i = point%y(i)
       
        nx = point%nx(i)
        ny = point%ny(i)


        tx = ny
        ty = -nx

        do j = 1, point%ypos_nbhs(i)

            k = point%ypos_conn(i,j)

            x_k = point%x(k)
            y_k = point%y(k)

            delx = x_k - x_i
            dely = y_k - y_i

            dels = delx*tx + dely*ty
            deln = delx*nx + dely*ny

            dist = dsqrt(dels*dels + deln*deln)
            weights = dist**power

            dels_weights = dels*weights
            deln_weights = deln*weights

            sum_delx_sqr = sum_delx_sqr + dels*dels_weights
            sum_dely_sqr = sum_dely_sqr + deln*deln_weights

            sum_delx_dely = sum_delx_dely + dels*deln_weights


            qtilde_i = point%q(:,i) - 0.5d0*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))

            call venkat_limiter(qtilde_i, phi_i, i)
            call venkat_limiter(qtilde_k, phi_k, k)

            qtilde_i = point%q(:,i) - 0.5d0*phi_i*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*phi_k*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))


            call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
            call flux_Gyp(G_i, nx, ny, u1, u2, rho, pr)

            call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
            call flux_Gyp(G_k, nx, ny, u1, u2, rho, pr)

            sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
            sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

        enddo

        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
        one_by_det = 1.d0/det

        G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det


    end subroutine




!	This subroutine evaluates the interior flux derivative dGx_neg


    subroutine interior_dGy_neg(G, i)


        implicit none

        integer :: i, j, k, r
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

        x_i = point%x(i)
        y_i = point%y(i)

        nx = point%nx(i)
        ny = point%ny(i)

        tx = ny
        ty = -nx

        do j = 1, point%yneg_nbhs(i)

            k = point%yneg_conn(i,j)

            x_k = point%x(k)
            y_k = point%y(k)

            delx = x_k - x_i
            dely = y_k - y_i

            dels = delx*tx + dely*ty    
            deln = delx*nx + dely*ny

            dist = dsqrt(dels*dels + deln*deln)
            weights = dist**power

            dels_weights = dels*weights
            deln_weights = deln*weights

            sum_delx_sqr = sum_delx_sqr + dels*dels_weights
            sum_dely_sqr = sum_dely_sqr + deln*deln_weights

            sum_delx_dely = sum_delx_dely + dels*deln_weights

            qtilde_i = point%q(:,i) - 0.5d0*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))

            call venkat_limiter(qtilde_i, phi_i, i)
            call venkat_limiter(qtilde_k, phi_k, k)

            qtilde_i = point%q(:,i) - 0.5d0*phi_i*(delx*point%dq(1,:,i) + dely*point%dq(2,:,i)) 
            qtilde_k = point%q(:,k) - 0.5d0*phi_k*(delx*point%dq(1,:,k) + dely*point%dq(2,:,k))

            call qtilde_to_primitive(qtilde_i, u1, u2, rho, pr)
            call flux_Gyn(G_i, nx, ny, u1, u2, rho, pr)

            call qtilde_to_primitive(qtilde_k, u1, u2, rho, pr)
            call flux_Gyn(G_k, nx, ny, u1, u2, rho, pr)

            sum_delx_delf = sum_delx_delf + (G_k - G_i)*dels_weights
            sum_dely_delf = sum_dely_delf + (G_k - G_i)*deln_weights

        enddo
    
        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
        one_by_det = 1.d0/det

        G = (sum_dely_delf*sum_delx_sqr - sum_delx_delf*sum_delx_dely)*one_by_det


    end subroutine


end module interior_fluxes_mod
