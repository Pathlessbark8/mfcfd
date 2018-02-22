#Start of the makefile

# defining the variables
objects=meshfree_solver.o data_structure.o point_preprocessor.o parameter.o
f90comp=gfortran  
ccomp=gcc

switch = -O3
execname: $(objects) 
	$(f90comp) -o execname $(switch) $(objects)

meshfree_solver.o:  point_preprocessor.o data_structure.o parameter.o
	$(f90comp) -c $(switch) meshfree_solver.f90

point_preprocessor.o: point_preprocessor.f90 data_structure.o
	$(f90comp) -c $(switch) point_preprocessor.f90

data_structure.o: data_structure.f90 parameter.o
	$(f90comp) -c $(switch) data_structure.f90

parameter.o: parameter.f90 
	$(f90comp) -c $(switch) parameter.f90


# Cleaning 

clean:
	rm -rf *.mod execname *.o

# End of makefile ..


