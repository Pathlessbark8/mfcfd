module compute_delta_mod

      use device_data_structure_mod

contains

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

end module
