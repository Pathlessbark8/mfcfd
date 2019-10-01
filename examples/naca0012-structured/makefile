ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_adjoint

all: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp adjoint/* . && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built mpi adjoint solver"

input:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile 
	cp mfpre/graphpart .
	rm -rf mfpre

clean:
	rm -rf src
	rm -f execname
