#Start of the makefile
include makefile.in

# defining the variables
objects=meshfree_solver.o data_structure.o point_preprocessor.o parameter.o initial_conditions.o q_lskum.o point_normals.o test.o generate_connectivity.o fpi_solver.o q_variables.o 
CC=gcc


switch = -O3

execname: $(objects) 
	$(FC)  -o execname $(switch) $(objects)

meshfree_solver.o:  point_preprocessor.o data_structure.o parameter.o initial_conditions.o q_lskum.o
	$(FC) -c $(switch) meshfree_solver.f90 

point_preprocessor.o: point_preprocessor.f90 data_structure.o
	$(FC) -c $(switch) point_preprocessor.f90

data_structure.o: data_structure.f90 parameter.o
	$(FC) -c $(switch) data_structure.f90

parameter.o: parameter.f90 
	$(FC) -c $(switch) parameter.f90

initial_conditions.o: initial_conditions.f90 data_structure.o 
	$(FC) -c $(switch) initial_conditions.f90

q_lskum.o: q_lskum.f90 data_structure.o point_normals.o  generate_connectivity.o fpi_solver.o 
	$(FC) -c $(switch) q_lskum.f90

point_normals.o: point_normals.f90 data_structure.o  
	$(FC) -c $(switch) point_normals.f90

test.o: test.f90 
	$(FC) -c $(switch) test.f90

generate_connectivity.o: generate_connectivity.f90 data_structure.o
	$(FC) -c $(switch) generate_connectivity.f90

fpi_solver.o: fpi_solver.f90 data_structure.o q_variables.o
	$(FC) -c $(switch) fpi_solver.f90

q_variables.o: q_variables.f90 data_structure.o
	$(FC) -c $(switch) q_variables.f90



# Cleaning 

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


