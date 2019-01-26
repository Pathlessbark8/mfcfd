program meshfree_solver

        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use initial_conditions_mod
        use q_lskum_mod
        use compute_force_coeffs_mod


        implicit none
        real*8  :: start,finish
        real*8  :: startr,finishr
        
        call cpu_time(start)

!       Set up case input

        call readinp()

!	Reading the input data ..

        call read_input_point_data()

!       Allocate solution variables

        call allocate_soln()

        plen = max_points
        print*,'No of points:',plen

!	Assign the initial conditions for the primitive variables ..	

        call initial_conditions()
        print*,'solution initialised'
       
!	Primal fixed point iterative solver ..
        
        call cpu_time(startr) 
        call q_lskum()
        call cpu_time(finishr) 

!       Save solution one last time
        call print_primal_output()


!       Deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()

        call cpu_time(finish) 
        print*,'Simulation completed!!!'
        print*,'runtime:',finishr-startr
        print*,'totaltime:',finish-start
        
end program meshfree_solver
