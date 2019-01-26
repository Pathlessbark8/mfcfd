#Start of the makefile

include makefile.in

# defining the variables

objects=meshfree_solver.o parameter.o data_structure.o point_preprocessor.o initial_conditions.o q_lskum.o point_normals.o  generate_connectivity.o device_data_structure.o post_processing.o q_variables.o compute_delta.o flux_residual.o state_update.o

execname: $(objects)  
	$(FC)  -o execname $(FFLAGS) $(objects)


meshfree_solver.o: meshfree_solver.F90 point_preprocessor.o data_structure.o parameter.o initial_conditions.o q_lskum.o post_processing.o
	$(FC) $(FFLAGS) -c  meshfree_solver.F90 

parameter.o: parameter.F90
	$(FC) $(FFLAGS)  -c  parameter.F90

point_preprocessor.o: point_preprocessor.F90 data_structure.o device_data_structure.o
	$(FC) $(FFLAGS) -c  point_preprocessor.F90 

data_structure.o: data_structure.F90 parameter.o
	$(FC) $(FFLAGS) -c  data_structure.F90

device_data_structure.o: device_data_structure.F90 parameter.o data_structure.o
	$(FC) $(FFLAGS) -c  device_data_structure.F90

initial_conditions.o: initial_conditions.F90 data_structure.o parameter.o
	$(FC) $(FFLAGS) -c  initial_conditions.F90

q_lskum.o: q_lskum.F90 data_structure.o point_normals.o generate_connectivity.o q_variables.o compute_delta.o state_update.o
	$(FC) $(FFLAGS) -c  q_lskum.F90

point_normals.o: point_normals.F90 data_structure.o
	$(FC) $(FFLAGS) -c  point_normals.F90

generate_connectivity.o: generate_connectivity.F90 data_structure.o 
	$(FC) $(FFLAGS) -c  generate_connectivity.F90

post_processing.o: post_processing.F90 data_structure.o
	$(FC) $(FFLAGS) -c  post_processing.F90

q_variables.o: q_variables.F90 device_data_structure.o
	$(FC) $(FFLAGS) -c  q_variables.F90

compute_delta.o: compute_delta.F90 device_data_structure.o
	$(FC) $(FFLAGS) -c  compute_delta.F90

flux_residual.o: flux_residual.F90 device_data_structure.o
	$(FC) $(FFLAGS) -c  flux_residual.F90

state_update.o: state_update.F90 device_data_structure.o
	$(FC) $(FFLAGS) -c  state_update.F90

clean::
	rm -rf *.mod execname *.o
