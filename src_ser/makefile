#Start of the makefile
include makefile.in

# defining the variables
objects=meshfree_solver.o data_structure.o point_preprocessor.o parameter.o initial_conditions.o
CC=gcc


switch = -O3

execname: $(objects) 
	$(FC)  -o execname $(switch) $(objects)

meshfree_solver.o:  point_preprocessor.o data_structure.o parameter.o initial_conditions.o
	$(FC) -c $(switch) meshfree_solver.f90 

point_preprocessor.o: point_preprocessor.f90 data_structure.o
	$(FC) -c $(switch) point_preprocessor.f90

data_structure.o: data_structure.f90 parameter.o
	$(FC) -c $(switch) data_structure.f90

parameter.o: parameter.f90 
	$(FC) -c $(switch) parameter.f90

initial_conditions.o: initial_conditions.f90 data_structure.o 
	$(FC) -c $(switch) initial_conditions.f90

# Cleaning 

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


