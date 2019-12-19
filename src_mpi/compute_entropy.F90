module compute_entropy_mod

        use data_structure_mod

        contains


                subroutine compute_entropy()


                        implicit none

                        integer :: k
			real*8 :: temp1, temp2
                        real*8 :: gtotal_entropy

                        total_entropy = 0.d0

                        temp2 = dlog(pr_inf)

                        do k = 1, local_points
                                temp1 = (point%prim(1,k))**gamma
                                temp1 = point%prim(4,k)/temp1
                                temp1 = dlog(temp1)
                                point%entropy(k) = dabs(temp1 - temp2)
                                
                                total_entropy = total_entropy + dabs(temp1 - temp2)
                        enddo

                        !call MPI_Reduce(total_entropy, gtotal_entropy , 1, &
                        !        & MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)

                        if(rank == 0) then
                                !write(*,*)"total entropy :", gtotal_entropy
                        end if

                end subroutine 

end module compute_entropy_mod
