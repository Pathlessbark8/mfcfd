ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_mpi_serial

all:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && make -f makefile.primal clean
	cp parameter.F90 src/
	cd src && make -j -f makefile.primal execname
	mv src/execname .
	@echo "Built mpi solver"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src
	rm -f execname graphpart
