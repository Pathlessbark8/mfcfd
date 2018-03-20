#Start of the makefile

include makefile.in




# defining the variables
objects=meshfree_solver.o data_structure.o point_preprocessor.o parameter.o utilities.o initial_conditions.o q_lskum.o point_normals.o test.o generate_connectivity.o



execname: $(objects)  
	$(FC)  -o execname $(switch) $(objects) ${PETSC_LIB}

meshfree_solver.o: point_preprocessor.o data_structure.o parameter.o initial_conditions.o q_lskum.o
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

q_lskum.o: q_lskum.F90 data_structure.o point_normals.o generate_connectivity.o
	$(FC) -c $(switch) q_lskum.F90

point_normals.o: point_normals.F90 data_structure.o
	$(FC) -c $(switch) point_normals.F90

test.o: test.F90 data_structure.o 
	$(FC) -c $(switch) test.F90

generate_connectivity.o: generate_connectivity.F90 data_structure.o 
	$(FC) -c $(switch) generate_connectivity.F90



# Cleaning 

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


