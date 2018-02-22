program main
    
    use data_structure_mod
	use point_preprocessor_mod

    implicit none

	call read_input_point_data()

    print*,point(:)%global_id

end program main