module q_lskum_mod

        use data_structure_mod
        use point_normals_mod    
        use generate_connectivity_mod
        use fpi_solver_mod
        use initial_conditions_mod

contains

        subroutine q_lskum()

                implicit none

                integer :: i

                call initial_conditions()

                call compute_normals()
                call generate_connectivity()

                ! Set U_old to U for first iteration
                do i=1,local_points
                        point%U_old(1,i) = point%prim(1,i)
                        point%U_old(2,i) = point%prim(1,i)*point%prim(2,i)
                        point%U_old(3,i) = point%prim(1,i)*point%prim(3,i)
                        point%U_old(4,i) = 2.5d0*point%prim(4,i) + 0.5d0*point%prim(1,i)*&
                                &(point%prim(2,i)*point%prim(2,i) +&
                                &point%prim(3,i)*point%prim(3,i))
                end do

                if(restart == 0)itr = 0
                
                do it = itr+1, itr+max_iters                
                        call fpi_solver(it)
                enddo
                
        end subroutine

end module q_lskum_mod
