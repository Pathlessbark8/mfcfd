#Start of the makefile

include makefile.in




# defining the variables
objects=meshfree_solver.o data_structure.o point_preprocessor.o parameter.o utilities.o



execname: $(objects)  
	$(FC)  -o execname $(switch) $(objects) ${PETSC_LIB}

meshfree_solver.o: point_preprocessor.o data_structure.o parameter.o
	$(FLINKER) -c  meshfree_solver.F90 

point_preprocessor.o: point_preprocessor.F90 data_structure.o
	$(FLINKER)  -c  point_preprocessor.F90 

data_structure.o: data_structure.F90 parameter.o 
	$(FC)  -c $(switch) data_structure.F90 

parameter.o: parameter.F90 
	$(FC)  -c $(switch) parameter.F90 

utilities.o: utilities.F90
	$(FC) -c $(switch) utilities.F90


# Cleaning 

clean::
	rm -rf *.mod execname *.o

# End of makefile ..


