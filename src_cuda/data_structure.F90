module data_structure_mod


        use parameter_mod

        implicit none

        integer :: wall_points,interior_points,outer_points

        type :: points

		real*8, dimension(:), allocatable :: x,y
                integer, dimension(:), allocatable :: local_id
                integer, dimension(:), allocatable :: global_id
                integer, dimension(:), allocatable :: left,right
                integer, dimension(:), allocatable :: flag_1 ! stores location of point
                integer, dimension(:), allocatable :: flag_2 ! stores shape point belongs to 
                integer, dimension(:), allocatable :: nbhs
                integer, dimension(:,:), allocatable :: conn
                
                integer, dimension(:), allocatable :: qtdepth

		real*8, dimension(:), allocatable :: nx, ny

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
        integer,allocatable,dimension(:) :: outer_points_index
        integer,allocatable,dimension(:) :: interior_points_index

        !iterations
        integer :: it

        real*8, allocatable, dimension(:)  :: Cl, Cd, Cm

        real*8 :: total_entropy

    contains

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
