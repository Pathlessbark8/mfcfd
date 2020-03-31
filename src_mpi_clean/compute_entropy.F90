module compute_entropy_mod
!#include <petsc/finclude/petscsys.h>

        use data_structure_mod
        use petsc_data_structure_mod

        contains


                subroutine compute_entropy()


                        implicit none

                        integer :: k
			real*8 :: temp1, temp2
                        real*8 :: gtotal_entropy
                        ! PetscErrorCode :: ierr

                        total_entropy = 0.d0

                        temp2 = dlog(pr_inf)

                        do k = 1, local_points
                                temp1 = (point%prim(1,k))**gamma
                                temp1 = point%prim(4,k)/temp1
                                temp1 = dlog(temp1)
                                point%entropy(k) = (temp1 - temp2)**2
                                
                                total_entropy = total_entropy + (temp1 - temp2)**2
                        enddo

                        call MPI_Reduce(total_entropy, gtotal_entropy , 1, &
                                & MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)

                        if(rank == 0) then
                                write(*,*)"total entropy :", gtotal_entropy
                        end if

                end subroutine
                
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
                    real*8 :: gtotal_enstrophy
                    PetscErrorCode :: ierr
                    
                    total_enstrophy = 0.d0
                    
                    do i = 1, local_points
                        
                        
                        x_i = point%x(i)
                        y_i = point%y(i)
                        
                        sum_delx_sqr = 0.d0
                        sum_dely_sqr = 0.d0
                        sum_delx_dely = 0.d0
                        
                        sum_delx_delu1 = 0.d0
                        sum_dely_delu1 = 0.d0
                        sum_delx_delu2 = 0.d0
                        sum_dely_delu2 = 0.d0
                        
                        do k = 1, point%nbhs(i)
                            
                            nbh = point%conn(i,k)
                            
                            x_k = point%x(nbh)
                            y_k = point%y(nbh)
                            
                            delx = x_k - x_i
                            dely = y_k - y_i
                            
                            dist = dsqrt(delx*delx + dely*dely)
                            weights = dist**power
                            
                            sum_delx_sqr = sum_delx_sqr + delx*delx*weights
                            sum_dely_sqr = sum_dely_sqr + dely*dely*weights
                            
                            sum_delx_dely = sum_delx_dely + delx*dely*weights
                            
                            sum_delx_delu1 = sum_delx_delu1 + weights*delx*(point%prim(2,nbh) - point%prim(2,i))
                            sum_delx_delu2 = sum_delx_delu2 + weights*delx*(point%prim(3,nbh) - point%prim(3,i))
                            sum_dely_delu1 = sum_dely_delu1 + weights*dely*(point%prim(2,nbh) - point%prim(2,i))
                            sum_dely_delu2 = sum_dely_delu2 + weights*dely*(point%prim(3,nbh) - point%prim(3,i))
                            
                        enddo
                        
                        det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely
                        one_by_det = 1.0d0/det
                        
                        du2_dx = (sum_delx_delu2*sum_dely_sqr - sum_dely_delu2*sum_delx_dely)*one_by_det
                        du1_dy = (sum_dely_delu1*sum_delx_sqr - sum_delx_delu1*sum_delx_dely)*one_by_det
                        
                        temp = du2_dx - du1_dy
                        
                        point%vorticity(i) = temp
                        
                        point%vorticity_sqr(i) = temp*temp
                        
                        total_enstrophy = total_enstrophy + point%vorticity_sqr(i)
                        
                    enddo
                    
                    call MPI_Reduce(total_enstrophy, gtotal_enstrophy , 1, &
                    & MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
                    
                    if(rank == 0) then
                        write(*,*)"total enstrophy :", gtotal_enstrophy
                    end if
                    
                end subroutine

end module compute_entropy_mod
