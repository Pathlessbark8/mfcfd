subroutine readparam()
#include <petsc/finclude/petscsys.h>
        use petscsys
        use data_structure_mod
        implicit none
        PetscErrorCode :: ierr
        PetscBool :: set
        character(len=64)  :: string
        character(len=20)  :: fmt1, fmt2

        if(rank==0) print*,'Reading parameters from parameter.in'


        cfl = 0.0d0 ! Default cfl number
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                                 '-cfl',cfl,set,ierr); CHKERRQ(ierr)

        max_iters = 10000 ! Default iterations
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-max_iters',max_iters,set,ierr); CHKERRQ(ierr)


        mach = 0.0d0 ! Default mach no
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-mach',mach,set,ierr); CHKERRQ(ierr)

        aoa = 0.0d0 ! Default aoa
        call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-aoa',aoa,set,ierr); CHKERRQ(ierr)
        !calculate theta
        theta = aoa*pi/180.d0

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
        
        solution_restart = 0 ! Default : no restart
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                            '-solution_restart',&
                            solution_restart,set,ierr); CHKERRQ(ierr)


        ! Print paramaters to screen
        fmt1 = "(5x,a10,i14)"
        fmt2 = "(5x,a10,e14.4)"
        write(*,fmt1) 'max_iters  =', max_iters
        write(*,fmt2) 'cfl   =', cfl
        write(*,fmt2) 'mach   =', mach
        write(*,fmt2) 'aoa  =', aoa
        write(*,fmt1) 'shapes      =', shapes


end subroutine 
