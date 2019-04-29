module q_lskum_mod

        use q_variables_mod
        use flux_residual_mod
        use state_update_mod
        use data_structure_mod
        use device_data_structure_mod
        use point_normals_mod    
        use generate_connectivity_mod
        use post_processing_mod
        use cudafor

        real*8 :: sum_res_sqr, residue
        real*8, device :: temp
        real*8, device, allocatable :: sum_res_sqr_d(:)

contains

        subroutine q_lskum()

                implicit none

                ! Grid and block dim
                type(dim3) :: grid , tBlock
                integer :: istat, stat
                integer :: i

                OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

                if(format_file == 2 .or. format_file == 3) call compute_normals()
                call generate_connectivity()
                write(*,*)'%%%%-Normals and connectivity generated-%%%'
                write(*,*)

                ! Transfer from host device
                call host_to_device()
                
                write(*,*)'%%%%%%%%%-Host to device performed-%%%%%%%%'
                write(*,*)

                allocate(sum_res_sqr_d(max_points))

                tBlock = dim3 (blockx, blocky, blockz)
                grid = dim3(ceiling(real(max_points)/ tBlock%x), 1, 1)
                
                write(*,*)'%%%%%%%%%%%%%%%-GPU size info-%%%%%%%%%%%%%'
                write(*,*) 'number of threads per block:',blockx*blocky*blockz
                write(*,*) 'grid dimension:',grid
                write(*,*) 'thread block dimension:',tBlock

                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-Iterations begin-%%%%%%%%%%%%'
                write(*,*)

                if(solution_restart == 0) itr = 0

                do it = itr+1, itr+max_iters

                        call eval_q_variables<<<grid, tBlock>>>(point_d%prim, point_d%q)

                        call eval_q_derivatives<<<grid, tBlock>>>(point_d%x, point_d%nbhs, &
                                & point_d%conn, point_d%q, point_d%dq)

                        call cal_flux_residual<<<grid, tBlock>>>(point_d%x, point_d%nx, &
                                & point_d%flag, point_d%nbhs, point_d%conn, &
                                & point_d%xpos_nbhs, point_d%xneg_nbhs, point_d%ypos_nbhs, &
                                & point_d%yneg_nbhs, point_d%xpos_conn, point_d%xneg_conn,&
                                & point_d%ypos_conn, point_d%yneg_conn, point_d%prim,  &
                                & point_d%q, point_d%dq, point_d%flux_res)

                        call state_update<<<grid, tBlock>>>(point_d%x, point_d%nx, point_d%flag, &
                                & point_d%nbhs, point_d%conn, point_d%prim, point_d%flux_res, sum_res_sqr_d)

                        istat = cudaGetLastError() 

                        if (istat .ne. 0) then
                                print*, cudaGetErrorString(istat) 
                                stop istat 
                        endif
                        
                        ! Residue norm evaluation
                        
                        temp = 0.0
                        !$cuf kernel do <<< *, * >>>
                        do i = 1, mp_d
                                temp = temp + sum_res_sqr_d(i)
                        end do

                        sum_res_sqr = temp

                        residue = dsqrt(sum_res_sqr)/max_points

                        if (it .le. 2 .and. solution_restart == 0) then
                                res_old = residue
                                residue = 0.0d0
                        else
                                residue = dlog10(residue/res_old)
                        end if

                        write(*,*) "iterations:", it, "residue:", residue
                        write(301,*) it, residue

                        if(mod(it,savesol) == 0) then
                                write(*,*)'%%%%%%%%%%%%%-Saving solution-%%%%%%%%%%%%%'
                                call device_to_host()
                                call print_primal_output()
                        end if
                enddo

                it = it - 1
                
                call device_to_host()
                
                write(*,*)
                write(*,*)'%%%%%%%%%-Device to host performed-%%%%%%%%'
                write(*,*)

                stat = cudaDeviceSynchronize()
                if (stat .ne. 0) then
                        print*, cudaGetErrorString(stat) 
                        stop stat 
                endif

                CLOSE(UNIT=301)

        end subroutine

end module q_lskum_mod
