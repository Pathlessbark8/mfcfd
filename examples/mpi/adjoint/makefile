ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_adjoint

shape: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp shape_sensitivity/* . && mv point_preprocessor_shape.F90 point_preprocessor.F90 && mv data_structure_shape.F90 data_structure_b.F90 && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built mpi adjoint solver for pure shape sensitivities"

aoa: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp aoa_sensitivity/* . && mv point_preprocessor_aoa.F90 point_preprocessor.F90 && mv data_structure_aoa.F90 data_structure_b.F90 && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built mpi adjoint solver for aoa based shape sensitivities"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile 
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src
	rm -f execname
