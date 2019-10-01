#Start of the makefile

include makefile.in

# defining the variables

objects=meshfree_solver.o parameter.o data_structure.o point_preprocessor.o initial_conditions.o point_normals.o  generate_connectivity.o post_processing.o objective_function.o compute_force_coeffs.o adaptation_sensors.o compute_entropy.o cuda.o


TARGET = execname

all:
	cat device_data_structure.cuf q_variables.cuf flux_residual.cuf state_update.cuf  q_lskum.cuf> cuda.cuf
	make -f makefile $(TARGET)


$(TARGET): $(objects) 
	$(FC) -o execname $(FFLAGS) $(objects)


meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure.o parameter.o initial_conditions.o cuda.o post_processing.o objective_function.o adaptation_sensors.o
	$(FC) $(FFLAGS) -c  meshfree_solver.F90 

parameter.o: parameter.F90
	$(FC) $(FFLAGS)  -c  parameter.F90

point_preprocessor.o: point_preprocessor.F90 data_structure.o cuda.o
	$(FC) $(FFLAGS) -c  point_preprocessor.F90 

data_structure.o: data_structure.F90 parameter.o
	$(FC) $(FFLAGS) -c  data_structure.F90

adaptation_sensors.o: adaptation_sensors.F90 data_structure.o
	$(FC) $(FFLAGS) -c  adaptation_sensors.F90


initial_conditions.o: initial_conditions.F90 data_structure.o parameter.o
	$(FC) $(FFLAGS) -c  initial_conditions.F90

cuda.o: cuda.cuf initial_conditions.o generate_connectivity.o point_normals.o post_processing.o data_structure.o parameter.o 
	$(FC) $(FFLAGS) -c  cuda.cuf

point_normals.o: point_normals.F90 data_structure.o
	$(FC) $(FFLAGS) -c  point_normals.F90

generate_connectivity.o: generate_connectivity.F90 data_structure.o 
	$(FC) $(FFLAGS) -c  generate_connectivity.F90

post_processing.o: post_processing.F90 data_structure.o
	$(FC) $(FFLAGS) -c  post_processing.F90

objective_function.o: objective_function.F90 data_structure.o compute_force_coeffs.o compute_entropy.o
	$(FC) $(FFLAGS) -c  objective_function.F90

compute_force_coeffs.o: compute_force_coeffs.F90 data_structure.o 
	$(FC) $(FFLAGS) -c  compute_force_coeffs.F90

compute_entropy.o: compute_entropy.F90 data_structure.o 
	$(FC) $(FFLAGS) -c  compute_entropy.F90

clean::
	rm -rf *.mod execname *.o *.ptx poly cuda.cuf
