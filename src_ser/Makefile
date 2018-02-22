#Start of the makefile

# defining the variables
objects=main.o data_structure.o point_preprocessor.o 
f90comp=gfortran  
ccomp=gcc

switch = -O3
execname: $(objects) 
	$(f90comp) -o execname $(switch) $(objects)

main.o:  point_preprocessor.o data_structure.o
	$(f90comp) -c $(switch) main.f90

point_preprocessor.o: point_preprocessor.f90 data_structure.o
	$(f90comp) -c $(switch) point_preprocessor.f90

data_structure.o: data_structure.f90 
	$(f90comp) -c $(switch) data_structure.f90


# Cleaning 

clean:
	rm $(objects) *.mod execname *.o

# End of makefile ..


