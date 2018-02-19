#Start of the makefile

# defining the variables

objects=meshfree_solver.o q_lskum.o parameter.o data_structure.o point_preprocessor.o generate_connectivity.o initial_conditions.o \
	quadrant_fluxes.o wall_fluxes.o split_fluxes.o interior_fluxes.o outer_fluxes.o flux_residual.o \
	fpi_solver.o state_update.o miscellaneous.o q_variables.o objective_function.o \
	compute_force_coeffs.o compute_entropy.o compute_enstrophy.o \
	post_processing.o limiters.o point_normals.o  
#
#
#
f90comp=gfortran  
ccomp=gcc

switch = -O3
#switch = -O3 -unroll=4 -depend

#Makefile

execname: $(objects) 

	$(f90comp) -o execname $(switch) $(objects)

meshfree_solver.o: meshfree_solver.f90 parameter.o point_preprocessor.o initial_conditions.o q_lskum.o \
			 data_structure.o post_processing.o 
	$(f90comp) -c $(switch) meshfree_solver.f90

q_lskum.o: q_lskum.f90 data_structure.o point_normals.o generate_connectivity.o fpi_solver.o
	$(f90comp) -c $(switch) q_lskum.f90

point_preprocessor.o: point_preprocessor.f90 data_structure.o
	$(f90comp) -c $(switch) point_preprocessor.f90

point_normals.o : point_normals.f90 data_structure.o
	$(f90comp) -c $(switch) point_normals.f90

generate_connectivity.o: generate_connectivity.f90 data_structure.o
	$(f90comp) -c $(switch) generate_connectivity.f90

initial_conditions.o: initial_conditions.f90 data_structure.o
	$(f90comp) -c $(switch) initial_conditions.f90

data_structure.o: data_structure.f90 parameter.o
	$(f90comp) -c $(switch) data_structure.f90

parameter.o: parameter.f90 
	$(f90comp) -c $(switch) parameter.f90

#
#
#
#
quadrant_fluxes.o: quadrant_fluxes.f90 parameter.o 
	$(f90comp) -c $(switch) quadrant_fluxes.f90
#
wall_fluxes.o: wall_fluxes.f90 quadrant_fluxes.o split_fluxes.o data_structure.o q_variables.o limiters.o
	$(f90comp) -c $(switch) wall_fluxes.f90
#
split_fluxes.o: split_fluxes.f90 parameter.o 
	$(f90comp) -c $(switch) split_fluxes.f90
#
interior_fluxes.o: interior_fluxes.f90 split_fluxes.o data_structure.o q_variables.o limiters.o
	$(f90comp) -c $(switch) interior_fluxes.f90
#
outer_fluxes.o: outer_fluxes.f90 split_fluxes.o quadrant_fluxes.o q_variables.o limiters.o data_structure.o
	$(f90comp) -c $(switch) outer_fluxes.f90
#
flux_residual.o: flux_residual.f90 interior_fluxes.o wall_fluxes.o outer_fluxes.o data_structure.o parameter.o
	$(f90comp) -c $(switch) flux_residual.f90
#
fpi_solver.o: fpi_solver.f90 q_variables.o flux_residual.o data_structure.o state_update.o objective_function.o 
	$(f90comp) -c $(switch) fpi_solver.f90
#
state_update.o: state_update.f90 flux_residual.o data_structure.o
	$(f90comp) -c $(switch) state_update.f90
#
miscellaneous.o: miscellaneous.f90 
	$(f90comp) -c $(switch) miscellaneous.f90
#
objective_function.o: objective_function.f90 data_structure.o compute_force_coeffs.o compute_entropy.o compute_enstrophy.o
	$(f90comp) -c $(switch) objective_function.f90
#
compute_force_coeffs.o: compute_force_coeffs.f90 data_structure.o
	$(f90comp) -c $(switch) compute_force_coeffs.f90
#
compute_entropy.o: compute_entropy.f90 data_structure.o
	$(f90comp) -c $(switch) compute_entropy.f90
#
compute_enstrophy.o: compute_enstrophy.f90 data_structure.o
	$(f90comp) -c $(switch) compute_enstrophy.f90
#
#	Files related to second order accuracy by q-variables ..
#
#
q_variables.o: q_variables.f90 data_structure.o 
	$(f90comp) -c $(switch) q_variables.f90
#
#
limiters.o: limiters.f90 data_structure.o
	$(f90comp) -c $(switch) limiters.f90 
#	
#
#	Files related to post-processing ..
#
#
post_processing.o: post_processing.f90 data_structure.o
	$(f90comp) -c $(switch) post_processing.f90
#
#
# Cleaning 

clean:
	rm $(objects) *.mod execname *.o

# End of makefile ..


