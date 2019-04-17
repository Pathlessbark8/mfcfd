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
	mv src/poly .
	@echo "Built cpu solver"

gpu: 
	rm -rf src
	cp -r $(SRCDIRC) src
	cd src && make -f makefile clean
	cp parameter* limiter.h src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built gpu solver"

cp:
	python3 cpplot.py

tec:
	./poly preprocessorfile.poly
	rm -rf preprocessorfile.1.node preprocessorfile.1.poly primal.1.ele
	gfortran tec_conversion.F90
	./a.out
	rm a.out


clean:
	rm -rf src cpplot_*
	rm -f execname
