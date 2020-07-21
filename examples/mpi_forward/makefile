ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_mpi_forward

phi:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp phi_files/* . && make -f makefile.forward clean
	cp parameter.F90 src/
	cd src && make -j -f makefile.forward execname
	mv src/execname .
	@echo "Built mpi forward solver w.r.t phi"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src
	rm -f execname graphpart
