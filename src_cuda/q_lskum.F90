module q_lskum_mod

        use q_variables_mod
        use compute_delta_mod
        use flux_residual_mod
        use state_update_mod
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

end module q_lskum_mod
