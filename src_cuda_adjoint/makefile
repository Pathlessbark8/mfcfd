#Start of the makefile
include makefile.in

#AD stuff
AD = adBuffer.o adStack.o revolve.o

# defining the variables

objects=meshfree_solver.o parameter.o data_structure_b.o point_preprocessor.o initial_conditions.o point_normals.o  generate_connectivity.o \
		post_processing.o objective_function.o compute_force_coeffs.o adaptation_sensors.o compute_entropy.o stagnation_values.o cuda.o \
		file_ops.o fpi_solver_b.o q_variables_b.o flux_residual_b.o wall_fluxes_b.o \
		interior_fluxes_b.o outer_fluxes_b.o limiters_b.o state_update_b.o objective_function_b.o \
		compute_force_coeffs_b.o compute_entropy_b.o split_fluxes_b.o quadrant_fluxes_b.o stagnation_values_b.o 

adBuffer.o : ADFirstAidKit/adBuffer.f
	$(FC) -c ADFirstAidKit/adBuffer.f

revolve.o : revolve.f
	$(FC) -c revolve.f

adStack.o : ADFirstAidKit/adStack.c
	$(CC) -c ADFirstAidKit/adStack.c

TARGET = execname

all:
	cat device_data_structure.cuf q_variables.cuf flux_residual.cuf state_update.cuf q_lskum.cuf > cuda.cuf
	make -f makefile $(TARGET)


$(TARGET): $(objects) $(AD)
	$(FC) -o execname $(FFLAGS) $(objects) $(AD)

SOURCEFILES = $(wildcard $(SRC)

meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure_b.o parameter.o initial_conditions.o cuda.o post_processing.o objective_function.o adaptation_sensors.o file_ops.o
	$(FC) $(FFLAGS) -c  meshfree_solver.F90 

parameter.o: parameter.F90
	$(FC) $(FFLAGS)  -c  parameter.F90

point_preprocessor.o: point_preprocessor.F90 data_structure_b.o cuda.o
	$(FC) $(FFLAGS) -c  point_preprocessor.F90 

data_structure_b.o: data_structure_b.F90 parameter.o
	$(FC) $(FFLAGS) -c  data_structure_b.F90

adaptation_sensors.o: adaptation_sensors.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c  adaptation_sensors.F90

initial_conditions.o: initial_conditions.F90 data_structure_b.o parameter.o
	$(FC) $(FFLAGS) -c  initial_conditions.F90

cuda.o: cuda.cuf initial_conditions.o generate_connectivity.o point_normals.o post_processing.o data_structure_b.o parameter.o stagnation_values.o file_ops.o fpi_solver_b.o
	$(FC) $(FFLAGS) -c  cuda.cuf

point_normals.o: point_normals.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c  point_normals.F90

generate_connectivity.o: generate_connectivity.F90 data_structure_b.o 
	$(FC) $(FFLAGS) -c  generate_connectivity.F90

post_processing.o: post_processing.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c  post_processing.F90

objective_function.o: objective_function.F90 data_structure_b.o compute_force_coeffs.o compute_entropy.o
	$(FC) $(FFLAGS) -c  objective_function.F90

compute_force_coeffs.o: compute_force_coeffs.F90 data_structure_b.o 
	$(FC) $(FFLAGS) -c  compute_force_coeffs.F90

compute_entropy.o: compute_entropy.F90 data_structure_b.o 
	$(FC) $(FFLAGS) -c  compute_entropy.F90

stagnation_values.o: stagnation_values.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c stagnation_values.F90

file_ops.o: file_ops.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c file_ops.F90

fpi_solver_b.o: fpi_solver_b.F90 data_structure_b.o q_variables_b.o flux_residual_b.o state_update_b.o objective_function_b.o stagnation_values_b.o
	$(FC) $(FFLAGS) -c fpi_solver_b.F90

q_variables_b.o: q_variables_b.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c  q_variables_b.F90

flux_residual_b.o: flux_residual_b.F90 data_structure_b.o parameter.o interior_fluxes_b.o wall_fluxes_b.o outer_fluxes_b.o
	$(FC) $(FFLAGS) -c  flux_residual_b.F90

state_update_b.o: state_update_b.F90 data_structure_b.o flux_residual_b.o
	$(FC) $(FFLAGS) -c  state_update_b.F90

interior_fluxes_b.o: interior_fluxes_b.F90 data_structure_b.o q_variables_b.o limiters_b.o split_fluxes_b.o #split_fluxes_b.o
	$(FC) $(FFLAGS) -c  interior_fluxes_b.F90

wall_fluxes_b.o: wall_fluxes_b.F90 quadrant_fluxes_b.o data_structure_b.o q_variables_b.o limiters_b.o split_fluxes_b.o  #split_fluxes_b.o # quadrant_fluxes_b.o
	$(FC) $(FFLAGS) -c  wall_fluxes_b.F90

outer_fluxes_b.o: outer_fluxes_b.F90 quadrant_fluxes_b.o q_variables_b.o limiters_b.o data_structure_b.o split_fluxes_b.o #split_fluxes_b.o #quadrant_fluxes_b.o
	$(FC) $(FFLAGS) -c  outer_fluxes_b.F90

limiters_b.o: limiters_b.F90 data_structure_b.o
	$(FC) $(FFLAGS) -c  limiters_b.F90

quadrant_fluxes_b.o: quadrant_fluxes_b.F90 parameter.o
	$(FC) $(FFLAGS) -c  quadrant_fluxes_b.F90
	
split_fluxes_b.o: split_fluxes_b.F90 parameter.o
	$(FC) $(FFLAGS) -c  split_fluxes_b.F90

compute_force_coeffs_b.o: compute_force_coeffs_b.F90 data_structure_b.o data_structure_b.o #petsc_bata_structure.o
	$(FC) $(FFLAGS) -c  compute_force_coeffs_b.F90

compute_entropy_b.o: compute_entropy_b.F90 data_structure_b.o #petsc_bata_structure.o
	$(FC) $(FFLAGS) -c  compute_entropy_b.F90

objective_function_b.o: objective_function_b.F90 data_structure_b.o compute_entropy_b.o compute_force_coeffs_b.o #compute_enstrophy_b.o
	$(FC) $(FFLAGS) -c  objective_function_b.F90

stagnation_values_b.o: stagnation_values_b.F90 data_structure_b.o data_structure_b.o
	$(FC) $(FFLAGS) -c  stagnation_values_b.F90

clean::
	rm -rf *.mod execname *.o *.ptx poly cuda.cuf
