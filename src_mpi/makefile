#Start of the makefile

include makefile.in




# defining the variables
objects=meshfree_solver.o q_lskum.o parameter.o data_structure.o point_preprocessor.o generate_connectivity.o \
initial_conditions.o quadrant_fluxes.o wall_fluxes.o split_fluxes.o interior_fluxes.o outer_fluxes.o flux_residual.o \
fpi_solver.o q_variables.o limiters.o point_normals.o test.o utilities.o state_update.o




execname: $(objects)  
	$(FC)  -o execname $(switch) $(objects) ${PETSC_LIB}

meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure.o parameter.o initial_conditions.o q_lskum.o
	$(FLINKER) -c  meshfree_solver.F90 

point_preprocessor.o: point_preprocessor.F90 data_structure.o
	$(FLINKER)  -c  point_preprocessor.F90 

data_structure.o: data_structure.F90 parameter.o 
	$(FLINKER)  -c  data_structure.F90 

parameter.o: parameter.F90 
	$(FC)  -c $(switch) parameter.F90 

utilities.o: utilities.F90
	$(FC) -c $(switch) utilities.F90

initial_conditions.o: initial_conditions.F90 data_structure.o
	$(FC) -c $(switch) initial_conditions.F90

q_lskum.o: q_lskum.F90 data_structure.o point_normals.o generate_connectivity.o fpi_solver.o
	$(FC) -c $(switch) q_lskum.F90

point_normals.o: point_normals.F90 data_structure.o
	$(FC) -c $(switch) point_normals.F90

generate_connectivity.o: generate_connectivity.F90 data_structure.o 
	$(FC) -c $(switch) generate_connectivity.F90

test.o: test.F90 data_structure.o 
	$(FC) -c $(switch) test.F90

fpi_solver.o: fpi_solver.F90 data_structure.o q_variables.o flux_residual.o state_update.o
	$(FLINKER) -c fpi_solver.F90

q_variables.o: q_variables.F90 data_structure.o
	$(FC) -c $(switch) q_variables.F90

flux_residual.o: flux_residual.F90 data_structure.o parameter.o interior_fluxes.o wall_fluxes.o outer_fluxes.o
	$(FC) -c $(switch) flux_residual.F90

quadrant_fluxes.o: quadrant_fluxes.F90 parameter.o
	$(FC) -c $(switch) quadrant_fluxes.F90

wall_fluxes.o: wall_fluxes.F90 quadrant_fluxes.o split_fluxes.o data_structure.o q_variables.o limiters.o
	$(FC) -c $(switch) wall_fluxes.F90

split_fluxes.o: split_fluxes.F90 parameter.o
	$(FC) -c $(switch) split_fluxes.F90

interior_fluxes.o: interior_fluxes.F90 split_fluxes.o data_structure.o q_variables.o limiters.o
	$(FC) -c $(switch) interior_fluxes.F90

outer_fluxes.o: outer_fluxes.F90 split_fluxes.o quadrant_fluxes.o q_variables.o limiters.o data_structure.o
	$(FC) -c $(switch) outer_fluxes.F90

state_update.o: state_update.F90 data_structure.o flux_residual.o
	$(FC) -c $(switch) state_update.F90
	
# Cleaning 

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


