module q_lskum_mod

        use q_variables_mod
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

                OPEN(UNIT=301,FILE="residue",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

                if(old_format == 0) call compute_normals()
                call generate_connectivity()

                ! Transfer from host device

                call host_to_device()

                tBlock = dim3 (blockx, blocky, blockz)
                grid = dim3(ceiling(real(max_points)/ tBlock%x), 1, 1)
                
                write(*,*)'%%%%%%%%%%%%%%%-GPU size info-%%%%%%%%%%%%%'
                write(*,*) 'number of threads per block:',blockx*blocky*blockz
                write(*,*) 'grid dimension:',grid
                write(*,*) 'thread block dimension:',tBlock

                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-Iterations begin-%%%%%%%%%%%%'
                write(*,*)

                do it = 0, max_iters
                        
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
                                & point_d%nbhs, point_d%conn, point_d%prim, point_d%flux_res)
                        
                        istat = cudaGetLastError() 
                        
                        if (istat .ne. 0) then
                                print*, cudaGetErrorString(istat) 
                                stop istat 
                        endif

                        write(*,'(a12,i8,a15,e30.20)')'iterations:',it
                        write(301, *) it

                enddo
                
                call device_to_host()

                CLOSE(UNIT=301)

        end subroutine

end module q_lskum_mod
