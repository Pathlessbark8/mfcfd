subroutine readcase()
! #include <petsc/finclude/petscsys.h>
    use petscsys
    use petsc_data_structure_mod
    use data_structure_mod
    implicit none
        ! PetscErrorCode :: ierr
    PetscBool :: set
    character(len=64) :: format_file, time, limiter, restart_solution
    character(len=64) :: solution_accuracy, tscheme, run_option

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

    tfinal = 1.0d20 ! Default final time: large value
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                '-tfinal',tfinal,set,ierr); CHKERRQ(ierr)

    time = 'steady' ! Default steady
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                  '-time',time,set,ierr); CHKERRQ(ierr)

    if(trim(time) == 'steady') then
        timestep = 0
    elseif(trim(time) == 'unsteady') then
        timestep = 1
    end if

    run_option = 'normal' ! Default normal run 
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                  '-run_option',run_option,set,ierr); CHKERRQ(ierr)

    if(trim(run_option) == 'normal') then
        runop = 1
    elseif(trim(run_option) == 'petsc') then
        runop = 2
    end if
    
    tscheme = 'ssprk43' ! Default first order in time
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                  '-tscheme',tscheme,set,ierr); CHKERRQ(ierr)

    if(trim(tscheme) == 'first') then
        rks = 1
        euler = 2.0d0
    elseif(trim(tscheme) == 'ssprk43') then
        rks = 4
        euler = 1.0d0
    end if

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

    vl_const = 150.0d0 ! Default VK limiter constant
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                '-vl_const',vl_const,set,ierr); CHKERRQ(ierr)
   
    restart_solution = 'no' ! Default : use initial conditions
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                  '-restart_solution',restart_solution,set,ierr); CHKERRQ(ierr)
   
    if(trim(restart_solution) == 'same') then
        restart = 1
    elseif(trim(restart_solution) == 'no') then
        restart = 0
    end if

    inner_iterations = 0
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
        '-inner_iterations',inner_iterations,set,ierr); CHKERRQ(ierr)

    interior_points_normal_flag = 0 ! Default : use nx as 0.0 and ny as 1.0
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                '-interior_points_normal_flag', &
                & interior_points_normal_flag,set,ierr); CHKERRQ(ierr)
   
    shapes = 1 ! Default : one shape
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                '-shapes',shapes,set,ierr); CHKERRQ(ierr)
    
    nsave = 1000000 ! Default : saving solution count
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                '-nsave', nsave,set,ierr); CHKERRQ(ierr)
    
    solution_accuracy = 'second' ! Default: second order
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,&
                  '-solution_accuracy',solution_accuracy,set,ierr); CHKERRQ(ierr)
    
    if(trim(solution_accuracy) == 'second') then
        fo_flag = 1.0d0
    elseif(trim(solution_accuracy) == 'first') then
        fo_flag = 0.0d0
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
        write(*,*) 'max_iters       :', max_iters
        write(*,*) 'cfl         :', cfl
        write(*,*) 'shapes          :', shapes
        write(*,*) 'timestep        :', time
        write(*,*) 'run_option      :', run_option
        write(*,*) 'tscheme         :', tscheme
        write(*,*) 'file format     :', format_file
        write(*,*) 'accuracy        :', solution_accuracy
        write(*,*) 'nsave           :', nsave
        write(*,*) 'power           :', power
        write(*,*) 'limiter_flag    :', limiter_flag
        write(*,*) 'inner_iterations    :', inner_iterations
        write(*,*) 'mach        :', mach
        write(*,*) 'aoa         :', aoa
        if(limiter_flag == 1)write(*,*) 'vl_const     :', vl_const
    end if

end subroutine 