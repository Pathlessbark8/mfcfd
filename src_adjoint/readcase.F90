subroutine readcase()
#include <petsc/finclude/petscsys.h>
        use petscsys
        use petsc_data_structure_mod
        use data_structure_mod_diff
        implicit none
        PetscErrorCode :: ierr
        PetscBool :: set
        character(len=64)  :: string

        cfl = 0.0d0 ! Default cfl number
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                                 '-cfl',cfl,set,ierr); CHKERRQ(ierr)

        max_iters = 10000000 ! Default iterations
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-max_iters',max_iters,set,ierr); CHKERRQ(ierr)

        power = 0.0d0 ! Default power
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-power',power,set,ierr); CHKERRQ(ierr)

        limiter_flag = 1 ! Default limiter => VK
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-limiter_flag',limiter_flag,set,ierr); CHKERRQ(ierr)

        vl_const = 150.0d0 ! Default VK limiter constant
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-vl_const',vl_const,set,ierr); CHKERRQ(ierr)
   
        initial_conditions_flag = 0 ! Default : use initial conditions
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-initial_conditions_flag',initial_conditions_flag,set,ierr)&
                            ; CHKERRQ(ierr)
   
        interior_points_normal_flag = 0 ! Default : use nx as 0.0 and ny as 1.0
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-interior_points_normal_flag',&
                            interior_points_normal_flag,set,ierr); CHKERRQ(ierr)
   
        shapes = 1 ! Default : one shape
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-shapes',shapes,set,ierr); CHKERRQ(ierr)
        
        restart = 0 ! Default : no restart
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-restart', restart,set,ierr); CHKERRQ(ierr)

        nsave = 1000 ! Default : saving solution count
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-nsave',&
                            nsave,set,ierr); CHKERRQ(ierr)
        
        fo_flag = 1.0d0 ! Default: second order
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-fo_flag',fo_flag,set,ierr); CHKERRQ(ierr)
        
        ad_mode = 0 ! Default adjoint mode => blackbox
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-ad_mode',ad_mode,set,ierr); CHKERRQ(ierr)
        
        chkpts = 50 ! Default checkpoints
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-chkpts',chkpts,set,ierr); CHKERRQ(ierr)
        
        obj_flag = 0 ! Default: no objective function
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-obj_flag',obj_flag,set,ierr); CHKERRQ(ierr)

        ! Print paramaters to screen
        if (rank==0) then
                write(*,*) 'max_iters :', max_iters
                write(*,*) 'cfl       :', cfl
                write(*,*) 'shapes    :', shapes
                if(ad_mode == 1) then
                        write(*,*) 'chkpts    :', chkpts
                end if
                write(*,*)
        end if

end subroutine 
