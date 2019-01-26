subroutine readcase()
#include <petsc/finclude/petscsys.h>
        use petscsys
        use data_structure_mod
        implicit none
        PetscErrorCode :: ierr
        PetscBool :: set
        character(len=64)  :: string
        character(len=20)  :: fmt1, fmt2

        if(rank==0) print*,'Reading parameters from case.in'


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

        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-tscheme',string,set,ierr); CHKERRQ(ierr)
        if(trim(string) == 'explicit')then
                tscheme = 0
                if (rank==0) print*,"explicit time scheme being used"
        else if(trim(string) == 'implicit')then
                tscheme = 1
                if (rank==0) print*,"implicit time scheme being used"
        end if

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
        
        solution_restart = 0 ! Default : no restart
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-solution_restart',&
                            solution_restart,set,ierr); CHKERRQ(ierr)

        nsave = 1000 ! Default : saving solution count
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-nsave',&
                            nsave,set,ierr); CHKERRQ(ierr)
        
        obj_flag = 0 ! Default : no objective function
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-obj_flag',obj_flag,set,ierr)

        old_format = 0 ! Default : new format
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-old_format',&
                            old_format,set,ierr); CHKERRQ(ierr)

        fo_flag = 1.0d0 ! Default: second order
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-fo_flag',fo_flag,set,ierr); CHKERRQ(ierr)
        
        if(obj_flag==0) then
                if(rank==0)print*,"no objective function chosen"
        elseif(obj_flag==1) then
                if(rank==0)print*,"Cl is chosen as the objective function"
                cl_flag = 1.0d0
        elseif(obj_flag==2) then
                if(rank==0)print*,"Cd is chosen as the objective function"
                cd_flag = 1.0d0
        elseif(obj_flag==3) then
                if(rank==0)print*,"Cm is chosen as the objective function"
                cm_flag = 1.0d0
        elseif(obj_flag==4) then
                if(rank==0)print*,"Cl/Cd is chosen as the objective function"
                cl_cd_flag = 1.0d0
        elseif(obj_flag==5) then
                if(rank==0)print*,"Total entropy is chosen as the objective function"
                ent_flag = 1.0d0
        elseif(obj_flag==6) then
                if(rank==0)print*,"Total enstrophy is chosen as the objective function"
                ens_flag = 1.0d0
        end if

        ! Print paramaters to screen
        fmt1 = "(5x,a11,i14)"
        fmt2 = "(5x,a11,e14.4)"
        if (rank==0) then
                write(*,fmt1) 'max_iters =', max_iters
                write(*,fmt2) 'cfl   =', cfl
                write(*,fmt2) 'mach   =', mach
                write(*,fmt2) 'aoa  =', aoa
                write(*,fmt1) 'shapes  =', shapes
        end if


end subroutine 
