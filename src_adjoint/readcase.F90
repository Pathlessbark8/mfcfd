subroutine readcase()
#include <petsc/finclude/petscsys.h>
        use petscsys
        use petsc_data_structure_mod
        use data_structure_mod_diff
        implicit none
        PetscErrorCode :: ierr
        PetscBool :: set
        character(len=64) :: format_file, time, limiter, restart_solution, tscheme
        character(len=64) :: solution_accuracy, adjoint_mode, objective_function

        cfl = 0.0d0 ! Default cfl number
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                                 '-cfl',cfl,set,ierr); CHKERRQ(ierr)

        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                                 '-mach',mach,set,ierr); CHKERRQ(ierr)
        
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                                 '-aoa',aoa,set,ierr); CHKERRQ(ierr)

        max_iters = 10000000 ! Default iterations
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-max_iters',max_iters,set,ierr); CHKERRQ(ierr)

        power = 0.0d0 ! Default power
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-power',power,set,ierr); CHKERRQ(ierr)
        limiter = 'venkat' ! Default limiter => VK
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-limiter',limiter,set,ierr); CHKERRQ(ierr)

        if(trim(limiter) == 'venkat') then
                limiter_flag = 1
        elseif(trim(limiter) == 'minmax') then
                limiter_flag = 2
        end if

        tscheme = 'first' ! Default first order in time
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-tscheme',tscheme,set,ierr); CHKERRQ(ierr)

        if(trim(tscheme) == 'first') then
                rks = 1
                euler = 2.0d0
        elseif(trim(time) == 'ssprk43') then
                rks = 4
                euler = 1.0d0
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

        restart_solution = 'no' ! Default : use initial conditions
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-restart_solution',restart_solution,set,ierr); CHKERRQ(ierr)

        if(trim(restart_solution) == 'same') then
                restart = 1
        elseif(trim(restart_solution) == 'no') then
                restart = 0
        end if
 
        nsave = 1000000 ! Default : saving solution count
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-nsave',&
                            nsave,set,ierr); CHKERRQ(ierr)

        solution_accuracy = 'second' ! Default: second order
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-solution_accuracy',solution_accuracy,set,ierr); CHKERRQ(ierr)

        if(trim(solution_accuracy) == 'second') then
                fo_flag = 1.0d0
        elseif(trim(solution_accuracy) == 'first') then
                fo_flag = 0.0d0
        end if

        adjoint_mode = 'checkpoints' ! Default : Checkpointing method
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-adjoint_mode',adjoint_mode,set,ierr); CHKERRQ(ierr)

        if(trim(adjoint_mode) == 'checkpoints') then
                ad_mode = 1
        elseif(trim(adjoint_mode) == 'blackbox') then
                ad_mode = 0
        end if
        
        chkpts = 50 ! Default checkpoints
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-chkpts',chkpts,set,ierr); CHKERRQ(ierr)
        
        objective_function = 'cl' ! Default: Cl
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-objective_function',objective_function, &
                            & set,ierr); CHKERRQ(ierr)

        if(trim(objective_function) == 'cl') then
                obj_flag = 1
        elseif(trim(objective_function) == 'cd') then
                obj_flag = 2
        elseif(trim(objective_function) == 'cm') then
                obj_flag = 3
        elseif(trim(objective_function) == 'totalentropy') then
                obj_flag = 4
        elseif(trim(objective_function) == 'totalenstrophy') then
                obj_flag = 5
        elseif(trim(objective_function) == 'clbycd') then
                obj_flag = 6
        elseif(trim(objective_function) == 'dcldalpha') then
                obj_flag = 7
        elseif(trim(objective_function) == 'dcddalpha') then
                obj_flag = 8
        elseif(trim(objective_function) == 'dcmdalpha') then
                obj_flag = 9
        else
                SETERRA(PETSC_COMM_WORLD,1,'Invalid objective function, check again')
        end if
        
        format_file = 'legacy' ! Default : legacy
        call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                              '-format_file',format_file,set,ierr); CHKERRQ(ierr)

        if(trim(format_file) == 'legacy') then
                format = 1
        elseif(trim(format_file) == 'quadtree') then
                format = 2
        end if

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
