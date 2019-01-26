ifndef FF_DIR
	$(error directory not set)
endif

SRCDIRM=$(FF_DIR)/src_mpi
SRCDIRC=$(FF_DIR)/src_cuda
SRCDIRS=$(FF_DIR)/src_ser

all:
	@echo "make cpu or gpu"

ser: 
	rm -rf src
	cp -r $(SRCDIRS) src
	cd src && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built serial solver"

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
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built gpu solver"

clean:
	rm -rf src
	rm -f execname
