ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_cuda_adjoint

phi:
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp phi_files/* . && make -f makefile clean
	cp parameter.F90 src/
	cd src && cat device_data_structure.cuf q_variables.cuf flux_residual.cuf state_update.cuf q_lskum.cuf > cuda.cuf
	cd src && make -j -f makefile execname
	mv src/execname .
	@echo "Built gpu adjoint solver w.r.t phi"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src point cpplot_*
	rm -f execname graphpart
