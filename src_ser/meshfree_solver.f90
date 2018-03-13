!
!	Updated on 26.12.2016
!
!
program meshfree_solver
!
!
	use parameter_mod
	use data_structure_mod
	use point_preprocessor_mod
    use initial_conditions_mod
!	use q_lskum_mod
!	use post_processing_mod
!
!
	implicit none
!	
!
!	Reading the input data ..
!
	call read_input_point_data()
!
!	Assign the initial conditions for the primitive variables ..	
!
	call initial_conditions()


        print*,'working till here'
!
!	Primal fixed point iterative solver ..
!
!	call q_lskum()
!
!	Printing the output (post-processing) ..
!
!	call print_primal_output()
!
!	
end program meshfree_solver
