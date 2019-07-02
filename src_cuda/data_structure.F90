module data_structure_mod


        use parameter_mod

        implicit none

        integer :: max_points
        integer :: wall_points,interior_points,outer_points,shape_points

        type :: points

		real*8, dimension(:), allocatable :: x,y
                integer, dimension(:), allocatable :: left,right
                integer, dimension(:), allocatable :: flag_1 ! stores location of point
                integer, dimension(:), allocatable :: flag_2 ! stores shape point belongs to 
                integer, dimension(:), allocatable :: nbhs
                integer, dimension(:,:), allocatable :: conn
                
                integer, dimension(:), allocatable :: qtdepth

		real*8, dimension(:), allocatable :: nx, ny

		real*8, dimension(:), allocatable :: min_dist

                real*8, dimension(:,:), allocatable :: prim

                real*8, dimension(:), allocatable :: sensor, D2_dist

                integer, dimension(:), allocatable :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer, dimension(:,:), allocatable :: xpos_conn, xneg_conn
                integer, dimension(:,:), allocatable :: ypos_conn, yneg_conn

                real*8, dimension(:), allocatable :: entropy
        end type points

 
        type(points) :: point

        save

        integer,allocatable,dimension(:) :: wall_points_index
        integer,allocatable,dimension(:) :: shape_points_index
        integer,allocatable,dimension(:) :: outer_points_index
        integer,allocatable,dimension(:) :: interior_points_index

        !iterations
        integer :: it, itr

        !Flag for time stepping
        integer :: rks = 1
        real*8 :: euler = 2.0d0
        character(len=20)  :: tscheme = 'first'

        real*8, allocatable, dimension(:)  :: Cl, Cd, Cm

        real*8 :: total_entropy

        real*8 :: res_old

        !The parameter CFL is the CFL number for stability ..
        real*8 :: CFL = 0.0d0

        integer :: max_iters = 10000000
!
!       The parameter power is used to specify the weights 
!       in the LS formula for the derivatives ..
!       power = 0.0d0, -2.0d0, -4.0d0, -6.0d0 ..
!       For example, power = -2.0 implies that
!       power = -2.0 => weights = 1/d^2
!       power = -4.0 => weights = 1/d^4
!
        real*8 :: power = 0.0d0
!
        real*8 :: VL_CONST = 150.0d0  ! Venkatakrishnan limiter constant ..

!       Restart solution parameter
        character(len=20)  :: restart_solution = 'no'
        integer :: solution_restart

!       format tag
        character(len=20)  :: format_file = 'legacy'
        integer :: file_format = 1

!       solution accuracy
        character(len=20)  :: solution_accuracy = 'second'
        real*8 :: f_o_flag

!       No of shapes
        integer :: shapes = 1

!       save frequency
        integer :: nsave = 10000000

!       Interior normal flag
        integer :: interior_points_normal_flag = 0

!       Block input
        integer :: blockx = 32, blocky = 1, blockz = 1

        namelist / input_parameters /   &
                      shapes, &
                             cfl, &
                       max_iters, &
                          blockx, &
                          blocky, &
                          blockz, &
                          vl_const, &
                          power, &
                          restart_solution, &
                          solution_accuracy, &
                          format_file, &
                          nsave, &
                          interior_points_normal_flag, &
                          tscheme, &
                          mach, &
                          aoa


    contains

        subroutine readnml()

                implicit none
                integer :: os


                open(unit=10,file='input.nml',form='formatted',status='old',iostat=os)
                read(unit=10,nml=input_parameters)

                close(10)
                write(*,*) 'Limiter:', 'venkat'
                write(*,*)
                write(*,*) '%%%%%%%%%%%%%%-Nml file info-%%%%%%%%%%%%%%'
                write(*,nml=input_parameters)
                write(*,*)

                if(solution_accuracy == 'second') then
                        f_o_flag = 1.0d0
                elseif(solution_accuracy == 'first') then
                        f_o_flag = 0.0d0
                end if

                if(restart_solution == 'no') then
                        solution_restart = 0
                elseif(restart_solution == 'yes') then
                        solution_restart = 1
                end if

                if(tscheme == 'first') then
                        rks = 1
                        euler = 2.0d0
                elseif(tscheme == 'ssprk43') then
                        rks = 4
                        euler = 1.0d0
                end if

                if(format_file == 'legacy') then
                        file_format = 1
                elseif(format_file == 'quadtree') then
                        file_format = 2
                        interior_points_normal_flag = 100000
                end if
        end subroutine


        subroutine allocate_soln()
                implicit none

                allocate(point%prim(4,max_points))

                allocate(point%sensor(max_points))
                allocate(point%D2_dist(max_points))
                
                allocate(point%xpos_nbhs(max_points))
                allocate(point%xneg_nbhs(max_points))
                allocate(point%ypos_nbhs(max_points))
                allocate(point%yneg_nbhs(max_points))

                allocate(point%xpos_conn(max_points,15))
                allocate(point%xneg_conn(max_points,15))

                allocate(point%ypos_conn(max_points,15))
                allocate(point%yneg_conn(max_points,15))

                allocate(Cl(shapes))
                allocate(Cd(shapes))
                allocate(Cm(shapes))

                allocate(point%entropy(max_points))
        end subroutine

        subroutine deallocate_soln()
                implicit none

                deallocate(point%prim)
                
                deallocate(point%sensor)
                deallocate(point%D2_dist)

                deallocate(point%xpos_nbhs)
                deallocate(point%xneg_nbhs)
                deallocate(point%ypos_nbhs)
                deallocate(point%yneg_nbhs)

                deallocate(point%xpos_conn)
                deallocate(point%xneg_conn)

                deallocate(point%ypos_conn)
                deallocate(point%yneg_conn)

                deallocate(Cl)
                deallocate(Cd)
                deallocate(Cm)

                deallocate(point%entropy)
        end subroutine

end module data_structure_mod
