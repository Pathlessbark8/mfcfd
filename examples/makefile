ifndef MF_DIR
	$(error directory not set)
endif

SRCDIRM=$(MF_DIR)/src_mpi
SRCDIRC=$(MF_DIR)/src_cuda

all:
	@echo "make cpu or gpu"

cpu: 
	rm -rf src
	cp -r $(SRCDIRM) src
	cd src && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built cpu solver"

gpu: 
	rm -rf src
	cp -r $(SRCDIRC) src
	cd src && make -f makefile clean
	cp parameter* limiter.h src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built gpu solver"

clean:
	rm -rf src
	rm -f execname
