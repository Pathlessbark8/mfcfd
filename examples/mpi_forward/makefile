# ifndef MF_DIR
	# $(error directory not set)
# endif

SRCDIR=/home/kumar/FORTRAN/mfcfd/src_mpi_forward

run:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp sensitivity_files/* . && make -f makefile.forward clean
	# cp parameter.F90 src/
	cd src && make -j -f makefile.forward execname
	mv src/execname .
	@echo "Built mpi forward solver w.r.t x y"

clean:
	rm -rf src
	rm -f execname graphpart
