module objective_function_mod
#include <petsc/finclude/petscsys.h>
        
        use petscsys
        use data_structure_mod
!	use compute_force_coeffs_mod
!	use compute_entropy_mod
        use compute_enstrophy_mod

        contains


                subroutine objective_function()


                        implicit none
                        real*8 :: lCFV
                        PetscErrorCode :: ierr

                        if(obj_func_flag==0) return


!                        call compute_cl_cd_cm()
                        !call compute_entropy()
                        call compute_enstrophy()
!
!			CFV = Cl_flag*Cl + Cd_flag*Cd + Cm_flag*Cm + Cl_Cd_flag*Cl/Cd
!			CFV = Cl/Cd
!			CFV = total_entropy
                        lCFV = total_enstrophy

                        call MPI_Reduce(lCFV, CFV , 1, MPI_DOUBLE, MPI_SUM, 0, &
                        PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

                        if(rank==0)write(*,'(a12,i8,a25,3e30.20)')'iterations:',it,'Objective function:',CFV



                end subroutine objective_function
!
!									
!
end module objective_function_mod
!		
!
!
