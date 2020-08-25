# ifndef MF_DIR
	# $(error directory not set)
# endif

SRCDIR=/home/kumar/FORTRAN/mfcfd/src_mpi_adjoint

phi:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp sensitivity_files/* . && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname
	mv src/execname .
	@echo "Built mpi adjoint solver w.r.t x y"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src
	rm -f execname graphpart