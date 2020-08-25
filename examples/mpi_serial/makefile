# ifndef MF_DIR
	# $(error directory not set)
# endif

SRCDIR=/home/kumar/FORTRAN/mfcfd/src_mpi_serial

all:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && make -f makefile.primal clean
	cp parameter.F90 src/
	cd src && make -j -f makefile.primal execname
	mv src/execname .
	@echo "Built mpi solver"

clean:
	rm -rf src
	rm -f execname graphpart
