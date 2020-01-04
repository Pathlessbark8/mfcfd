ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_forward

# shape:
# 	rm -rf src
# 	cp -r $(SRCDIR) src
# 	cd src && cp shape_sensitivity/* . && make -f makefile clean
# 	cp parameter.F90 src/
# 	cd src && make -j -f makefile execname
# 	mv src/execname .
# 	@echo "Built mpi forward solver w.r.t shape"

# aoa:
# 	rm -rf src
# 	cp -r $(SRCDIR) src
# 	cd src && cp aoa_sensitivity/* . && make -f makefile clean
# 	cp parameter.F90 src/
# 	cd src && make -j -f makefile execname
# 	mv src/execname .
# 	@echo "Built mpi forward solver w.r.t aoa"

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
	rm -rf src point
	rm -f execname graphpart
