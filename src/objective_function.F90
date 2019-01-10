module objective_function_mod
#include <petsc/finclude/petscsys.h>
        
        use petscsys
        use data_structure_mod
        use compute_force_coeffs_mod
        use compute_entropy_mod
        use compute_enstrophy_mod

        contains


                subroutine objective_function()


                        implicit none
                        real*8, dimension(shapes) :: lCFV
                        PetscErrorCode :: ierr

                        if(obj_flag.ne.0)OPEN(UNIT=401,FILE="cfv.dat")




                        call compute_cl_cd_cm()
                        call compute_entropy()
                        call compute_enstrophy()

                        lCFV = Cl_flag*Cl + Cd_flag*Cd + Cm_flag*Cm + Cl_Cd_flag*Cl/Cd &
                                + Ent_flag*total_entropy + Ens_flag*total_enstrophy

                        call MPI_Reduce(lCFV, CFV , shapes, MPI_DOUBLE, MPI_SUM, 0, &
                        PETSC_COMM_WORLD, ierr)

                        if(rank==0.and.obj_flag.ne.0) then
                                write(*,'(a12,i8,a25,3e30.20)')'iterations:',it,'Objective function:',CFV
                                write(401,*)it,CFV
                        end if



                end subroutine objective_function
!
!									
!
end module objective_function_mod
!		
!
!
