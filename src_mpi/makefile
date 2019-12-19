#Start of the makefile

include makefile.in




# defining the variables

objects=meshfree_solver.o parameter.o data_structure.o petsc_data_structure.o point_preprocessor.o utilities.o initial_conditions.o q_lskum.o point_normals.o  generate_connectivity.o fpi_solver.o q_variables.o flux_residual.o wall_fluxes.o interior_fluxes.o outer_fluxes.o split_fluxes.o limiters.o quadrant_fluxes.o state_update.o readcase.o post_processing.o objective_function.o compute_force_coeffs.o compute_entropy.o compute_enstrophy.o


execname: $(objects)  
	$(FC)  -o execname $(FFLAGS) $(objects) $(PETSC_LIB) $(LIBS)


meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure.o petsc_data_structure.o parameter.o q_lskum.o post_processing.o compute_force_coeffs.o
	$(FC) $(FFLAGS) -c  meshfree_solver.F90 

point_preprocessor.o: point_preprocessor.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  point_preprocessor.F90 

data_structure.o: data_structure.F90 parameter.o 
	$(FC) $(FFLAGS) -c  data_structure.F90 

petsc_data_structure.o: petsc_data_structure.F90 parameter.o data_structure.o
	$(FC) $(FFLAGS) -c  petsc_data_structure.F90 

parameter.o: parameter.F90 
	$(FC) $(FFLAGS)  -c  parameter.F90 

readcase.o: readcase.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  readcase.F90 

utilities.o: utilities.F90
	$(FC) $(FFLAGS) -c  utilities.F90

initial_conditions.o: initial_conditions.F90 data_structure.o parameter.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  initial_conditions.F90

q_lskum.o: q_lskum.F90 data_structure.o point_normals.o generate_connectivity.o fpi_solver.o initial_conditions.o
	$(FC) $(FFLAGS) -c  q_lskum.F90

point_normals.o: point_normals.F90 data_structure.o
	$(FC) $(FFLAGS) -c  point_normals.F90

generate_connectivity.o: generate_connectivity.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  generate_connectivity.F90

fpi_solver.o: fpi_solver.F90 data_structure.o q_variables.o flux_residual.o state_update.o objective_function.o post_processing.o 
	$(FC) $(FFLAGS) -c fpi_solver.F90

q_variables.o: q_variables.F90 data_structure.o
	$(FC) $(FFLAGS) -c  q_variables.F90

flux_residual.o: flux_residual.F90 data_structure.o parameter.o interior_fluxes.o wall_fluxes.o outer_fluxes.o
	$(FC) $(FFLAGS) -c  flux_residual.F90

quadrant_fluxes.o: quadrant_fluxes.F90 parameter.o
	$(FC) $(FFLAGS) -c  quadrant_fluxes.F90

wall_fluxes.o: wall_fluxes.F90 quadrant_fluxes.o split_fluxes.o data_structure.o q_variables.o limiters.o
	$(FC) $(FFLAGS) -c  wall_fluxes.F90

split_fluxes.o: split_fluxes.F90 parameter.o
	$(FC) $(FFLAGS) -c  split_fluxes.F90

interior_fluxes.o: interior_fluxes.F90 split_fluxes.o data_structure.o q_variables.o limiters.o
	$(FC) $(FFLAGS) -c  interior_fluxes.F90

outer_fluxes.o: outer_fluxes.F90 split_fluxes.o quadrant_fluxes.o q_variables.o limiters.o data_structure.o
	$(FC) $(FFLAGS) -c  outer_fluxes.F90

state_update.o: state_update.F90 data_structure.o flux_residual.o
	$(FC) $(FFLAGS) -c  state_update.F90

limiters.o: limiters.F90 data_structure.o
	$(FC) $(FFLAGS) -c  limiters.F90
	
objective_function.o: objective_function.F90 data_structure.o  compute_enstrophy.o compute_entropy.o compute_force_coeffs.o
	$(FC) $(FFLAGS) -c  objective_function.F90
	
compute_force_coeffs.o: compute_force_coeffs.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_force_coeffs.F90

compute_enstrophy.o: compute_enstrophy.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_enstrophy.F90

compute_entropy.o: compute_entropy.F90 data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  compute_entropy.F90

post_processing.o: data_structure.o petsc_data_structure.o
	$(FC) $(FFLAGS) -c  post_processing.F90

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


