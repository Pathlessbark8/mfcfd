module compute_enstrophy_mod
       
        use data_structure_mod

        contains


                subroutine compute_enstrophy()


                        implicit none
                        
                        
                        integer :: i, k, r, nbh
                
                        real*8 :: x_i, y_i, x_k, y_k
			real*8 :: delx, dely, dist, weights
			real*8 :: sum_delx_sqr, sum_dely_sqr, sum_delx_dely
			real*8 :: sum_delx_delu1, sum_delx_delu2, sum_dely_delu1, sum_dely_delu2
			real*8 :: det
			real*8 :: one_by_det
			real*8 :: du1_dy, du2_dx, temp


                        total_enstrophy = 0.d0

                        do i = 1, local_points
                        

                                x_i = p%x(i)
                                y_i = p%y(i)

                                sum_delx_sqr = 0.d0
                                sum_dely_sqr = 0.d0
                                sum_delx_dely = 0.d0

                                sum_delx_delu1 = 0.d0
                                sum_dely_delu1 = 0.d0
                                sum_delx_delu2 = 0.d0
                                sum_dely_delu2 = 0.d0

                                do k = 1, p%nbhs(i)

                                        nbh = p%conn(i,k)

                                        x_k = p%x(nbh)
                                        y_k = p%y(nbh)

                                        delx = x_k - x_i
                                        dely = y_k - y_i

                                        dist = dsqrt(delx*delx + dely*dely)
                                        weights = dist**power

                                        sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                                        sum_dely_sqr = sum_dely_sqr + dely*dely*weights

                                        sum_delx_dely = sum_delx_dely + delx*dely*weights

                                        sum_delx_delu1 = sum_delx_delu1 + weights*delx*(p%prim(2,nbh) - p%prim(2,i))
                                        sum_delx_delu2 = sum_delx_delu2 + weights*delx*(p%prim(3,nbh) - p%prim(3,i))
                                        sum_dely_delu1 = sum_dely_delu1 + weights*dely*(p%prim(2,nbh) - p%prim(2,i))
                                        sum_dely_delu2 = sum_dely_delu2 + weights*dely*(p%prim(3,nbh) - p%prim(3,i))

                                enddo

                                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                                one_by_det = 1.0d0/det

                                du2_dx = (sum_delx_delu2*sum_dely_sqr - sum_dely_delu2*sum_delx_dely)*one_by_det
                                du1_dy = (sum_dely_delu1*sum_delx_sqr - sum_delx_delu1*sum_delx_dely)*one_by_det

                                temp = du2_dx - du1_dy

                                p%vorticity(i) = temp

                                p%vorticity_sqr(i) = temp*temp

                                total_enstrophy = total_enstrophy + p%vorticity_sqr(i)


                        enddo


        end subroutine
end module compute_enstrophy_mod
