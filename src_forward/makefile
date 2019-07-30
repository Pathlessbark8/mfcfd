#Start of the makefile

include makefile.in




# defining the variables

objects=meshfree_solver.o parameter.o data_structure_d.o petsc_data_structure.o point_preprocessor.o utilities.o initial_conditions.o q_lskum_d.o point_normals_d.o  generate_connectivity_d.o fpi_solver_d.o q_variables_d.o flux_residual_d.o wall_fluxes_d.o interior_fluxes_d.o outer_fluxes_d.o split_fluxes_d.o limiters_d.o quadrant_fluxes_d.o state_update_d.o readcase.o objective_function_d.o compute_force_coeffs_d.o compute_entropy_d.o compute_enstrophy_d.o post_processing.o


execname: $(objects)  
	$(FC)  -o execname $(FFLAGS) $(objects) $(PETSC_LIB)


meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure_d.o petsc_data_structure.o parameter.o q_lskum_d.o post_processing.o
	$(FC) $(FFLAGS) -c  meshfree_solver.F90 

point_preprocessor.o: point_preprocessor.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  point_preprocessor.F90 

post_processing.o: post_processing.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  post_processing.F90 

data_structure_d.o: data_structure_d.F90 parameter.o 
	$(FC) $(FFLAGS) -c  data_structure_d.F90 

petsc_data_structure.o: petsc_data_structure.F90 parameter.o data_structure_d.o
	$(FC) $(FFLAGS) -c  petsc_data_structure.F90 

parameter.o: parameter.F90 
	$(FC) $(FFLAGS)  -c  parameter.F90 

readcase.o: readcase.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  readcase.F90 

utilities.o: utilities.F90
	$(FC) $(FFLAGS) -c  utilities.F90

initial_conditions.o: initial_conditions.F90 data_structure_d.o parameter.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  initial_conditions.F90

q_lskum_d.o: q_lskum_d.F90 data_structure_d.o point_normals_d.o generate_connectivity_d.o fpi_solver_d.o initial_conditions.o
	$(FC) $(FFLAGS) -c  q_lskum_d.F90

point_normals_d.o: point_normals_d.F90 data_structure_d.o
	$(FC) $(FFLAGS) -c  point_normals_d.F90

generate_connectivity_d.o: generate_connectivity_d.F90 data_structure_d.o 
	$(FC) $(FFLAGS) -c  generate_connectivity_d.F90

fpi_solver_d.o: fpi_solver_d.F90 data_structure_d.o q_variables_d.o flux_residual_d.o state_update_d.o objective_function_d.o  
	$(FC) $(FFLAGS) -c fpi_solver_d.F90

q_variables_d.o: q_variables_d.F90 data_structure_d.o
	$(FC) $(FFLAGS) -c  q_variables_d.F90

flux_residual_d.o: flux_residual_d.F90 data_structure_d.o parameter.o interior_fluxes_d.o wall_fluxes_d.o outer_fluxes_d.o
	$(FC) $(FFLAGS) -c  flux_residual_d.F90

quadrant_fluxes_d.o: quadrant_fluxes_d.F90 parameter.o
	$(FC) $(FFLAGS) -c  quadrant_fluxes_d.F90

wall_fluxes_d.o: wall_fluxes_d.F90 quadrant_fluxes_d.o split_fluxes_d.o data_structure_d.o q_variables_d.o limiters_d.o
	$(FC) $(FFLAGS) -c  wall_fluxes_d.F90

split_fluxes_d.o: split_fluxes_d.F90 parameter.o
	$(FC) $(FFLAGS) -c  split_fluxes_d.F90

interior_fluxes_d.o: interior_fluxes_d.F90 split_fluxes_d.o data_structure_d.o q_variables_d.o limiters_d.o
	$(FC) $(FFLAGS) -c  interior_fluxes_d.F90

outer_fluxes_d.o: outer_fluxes_d.F90 split_fluxes_d.o quadrant_fluxes_d.o q_variables_d.o limiters_d.o data_structure_d.o
	$(FC) $(FFLAGS) -c  outer_fluxes_d.F90

state_update_d.o: state_update_d.F90 data_structure_d.o flux_residual_d.o
	$(FC) $(FFLAGS) -c  state_update_d.F90

limiters_d.o: limiters_d.F90 data_structure_d.o
	$(FC) $(FFLAGS) -c  limiters_d.F90
	
objective_function_d.o: objective_function_d.F90 data_structure_d.o  compute_enstrophy_d.o compute_entropy_d.o compute_force_coeffs_d.o
	$(FC) $(FFLAGS) -c  objective_function_d.F90
	
compute_force_coeffs_d.o: compute_force_coeffs_d.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_force_coeffs_d.F90

compute_enstrophy_d.o: compute_enstrophy_d.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_enstrophy_d.F90

compute_entropy_d.o: compute_entropy_d.F90 data_structure_d.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_entropy_d.F90

clean::
	rm -rf *.mod execname *.o
