module file_ops_mod

    use data_structure_mod_diff

    namelist / input_parameters /   &
    shapes, &
        cfl, &
        chkpts, &
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
        adjoint_mode, &
        read_phi, &
        ! mach, &
        ! aoa, &
        inner_iterations

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

        if(adjoint_mode == 'checkpoints') then
            ad_mode = 1
        elseif(adjoint_mode == 'blackbox') then
            ad_mode = 0
        end if

        if(restart_solution == 'no') then
            restart = 0
        elseif(restart_solution == 'yes') then
            restart = 1
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

        if(read_phi == 'yes') then
            read_phi_file = 1
        elseif(read_phi == 'no') then
            read_phi_file = 0
        endif

    end subroutine

end module file_ops_mod