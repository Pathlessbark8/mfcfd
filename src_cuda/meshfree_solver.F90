program meshfree_solver

        use cudafor
        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use initial_conditions_mod
        use q_lskum_mod
        use post_processing_mod
        use objective_function_mod
        use adaptation_sensors_mod


        implicit none
        integer :: istat, i, nDevices=0
        real*8  :: start,finish, runtime
        type(cudaDeviceProp) :: prop
        
        call cpu_time(start)

        write(*,*)
        write(*,*)'%%%%%%%%-CUDA Fortran Meshfree Code-%%%%%%%'
        write(*,*)
        write(*,*)'%%%%%%%%%%%%%%%-Device info-%%%%%%%%%%%%%%%'
        istat = cudaGetDeviceCount ( nDevices )
        do i = 0, nDevices - 1
                istat = cudaGetDeviceProperties(prop, i)
                write(*,*)'Device Name:               ', trim(prop%name)
                write(*,*)'Compute Capability:        ',prop%major, prop%minor
                write(*,*)'Device number:             ',i
                write(*,*)'MemoryClockRate(KHz):      ',prop%memoryBusWidth
                write(*,*)'PeakMemoryBandwidth(GB/s): ',2.0 *prop%memoryClockRate * &
                        & (prop%memoryBusWidth/8) * 1.e-6
                write(*,*)
        end do

!       Set up case input
        call readnml()

!	Reading the input data ..
        write(*,*) '%%%%%%%%%%%%-Reading point file-%%%%%%%%%%%'
        call read_input_point_data()
        write(*,*) 'Number of points:         ', max_points
        write(*,*) 'Number of wall points:    ', wall_points
        write(*,*) 'Number of shape points:   ', shape_points
        write(*,*) 'Number of interior points:', interior_points
        write(*,*) 'Number of outer points:   ', outer_points
        write(*,*)

!       Allocate solution variables
        call allocate_soln()

!       Allocate device solution variables
        call allocate_device_soln()

!	Assign the initial conditions for the primitive variables ..	
        call initial_conditions()
        write(*,*) '%%%%%%%%%%%-Solution initialized-%%%%%%%%%%'
        write(*,*)
       
!	Primal fixed point iterative solver ..
        call q_lskum(runtime)

!       Compute sensor values
        call compute_adapt_sensor()

!       Objective function computation
        call objective_function()

!       Save solution one last time
        call print_primal_output()

!       Deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()
        call deallocate_device_soln()

        call cpu_time(finish) 
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',finish-start,'seconds'

end program meshfree_solver