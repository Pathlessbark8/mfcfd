module data_structure_mod


        use parameter_mod

        implicit none

        integer, parameter :: max_points=9600
        integer, parameter :: wall_points=159,interior_points=9281,outer_points=160,shape_points=160
        !       No of shapes
        integer, parameter :: shapes = 1

        type :: points

		real*8, dimension(max_points) :: x,y
                integer, dimension(max_points) :: left,right
                integer, dimension(max_points) :: flag_1 ! stores location of point
                integer, dimension(max_points) :: flag_2 ! stores shape point belongs to
                integer, dimension(max_points) :: nbhs
                real*8, dimension(max_points) :: delta
                integer, dimension(max_points,15) :: conn

                integer, dimension(max_points) :: qtdepth

		real*8, dimension(max_points) :: nx, ny

		real*8, dimension(max_points) :: min_dist

                real*8, dimension(4,max_points) :: prim
                real*8, dimension(4,max_points) :: prim_old

                real*8, dimension(4,max_points) :: U

                real*8, dimension(4,max_points) :: q
                real*8, dimension(4,max_points) :: flux_res

                real*8, dimension(2,4,max_points) :: qm
                real*8, dimension(2,4,max_points) :: dq
                real*8, dimension(3,4,max_points) :: ddq
                real*8, dimension(3,4,max_points) :: temp

                real*8, dimension(4,max_points) :: phi1, phi2

                real*8, dimension(max_points) :: sensor, D2_dist

                integer, dimension(max_points) :: xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs
                integer, dimension(max_points,15) :: xpos_conn, xneg_conn
                integer, dimension(max_points,15) :: ypos_conn, yneg_conn

                real*8, dimension(max_points) :: entropy

        end type points


        type(points) :: point

        save

        integer,dimension(wall_points) :: wall_points_index
        integer,dimension(shape_points) :: shape_points_index
        integer,dimension(outer_points) :: outer_points_index
        integer,dimension(interior_points) :: interior_points_index

        real*8 :: cost_func

        !iterations
        integer :: it, itr

        !Flag for time stepping
        integer :: rks = 1
        !checkpoints
        integer :: chkpts
        real*8 :: euler = 2.0d0
        character(len=20)  :: tscheme = 'first'

        real*8, dimension(shapes) :: Cl, Cd, Cm

        real*8 :: total_entropy

        real*8 :: res_old, res_new, residue

        real*8 :: max_res

        integer :: max_res_point

        real*8 :: sum_res_sqr

        !The parameter CFL is the CFL number for stability ..
        real*8 :: cfl = 0.0d0

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

!       Inner Iterations Loop count
        integer :: inner_iterations = 0

!       format tag
        character(len=20)  :: format_file = 'legacy'
        integer :: file_format = 1

!       solution accuracy
        character(len=20)  :: solution_accuracy = 'second'
        real*8 :: f_o_flag

!       save frequency
        integer :: nsave = 10000000

!       Interior normal flag
        integer :: interior_points_normal_flag = 0

!       Restart
        integer :: restart = 0

!       Block input
        integer :: blockx = 32, blocky = 1, blockz = 1

        ! contains

        ! subroutine allocate_soln()
        !         implicit none

        !         allocate(point%sensor(max_points))
        !         allocate(point%D2_dist(max_points))
        !         allocate(point%delta(max_points))

        !         allocate(point%xpos_nbhs(max_points))
        !         allocate(point%xneg_nbhs(max_points))
        !         allocate(point%ypos_nbhs(max_points))
        !         allocate(point%yneg_nbhs(max_points))

        !         allocate(point%xpos_conn(max_points,15))
        !         allocate(point%xneg_conn(max_points,15))

        !         allocate(point%ypos_conn(max_points,15))
        !         allocate(point%yneg_conn(max_points,15))


        !         ! allocate(point%U_old(4,max_points))
        !         allocate(point%U(4,max_points))

        !         allocate(point%prim(4,max_points))
        !         allocate(point%prim_old(4,max_points))

        !         allocate(point%flux_res(4,max_points))

        !         allocate(point%q(4,max_points))
        !         allocate(point%phi1(4,max_points))
        !         allocate(point%phi2(4,max_points))
        !         allocate(point%dq(2,4,max_points))
        !         allocate(point%ddq(3,4,max_points))
        !         allocate(point%temp(3,4,max_points))
        !         allocate(point%qm(2,4,max_points))

        !         allocate(Cl(shapes))
        !         allocate(Cd(shapes))
        !         allocate(Cm(shapes))

        !         allocate(point%entropy(max_points))
        ! end subroutine

        ! subroutine deallocate_soln()
        !         implicit none

        !         deallocate(point%prim)

        !         deallocate(point%sensor)
        !         deallocate(point%D2_dist)

        !         deallocate(point%xpos_nbhs)
        !         deallocate(point%xneg_nbhs)
        !         deallocate(point%ypos_nbhs)
        !         deallocate(point%yneg_nbhs)

        !         deallocate(point%xpos_conn)
        !         deallocate(point%xneg_conn)

        !         deallocate(point%ypos_conn)
        !         deallocate(point%yneg_conn)

        !         deallocate(Cl)
        !         deallocate(Cd)
        !         deallocate(Cm)

        !         deallocate(point%entropy)
        ! end subroutine

end module data_structure_mod