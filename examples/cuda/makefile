ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_cuda

gpu: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built gpu solver"

point:
	cp -r $(MF_DIR)/mfpre .
	cd mfpre && make -f makefile 
	cp mfpre/graphpart .
	rm -rf mfpre

cp:
	python3 cpplot.py

tec:
	cp src/poly .
	cp $(MF_DIR)/examples/cuda/tec_conversion.F90 .
	./poly preprocessorfile.poly
	rm -rf preprocessorfile.1.node preprocessorfile.1.poly primal.1.ele
	gfortran tec_conversion.F90
	./a.out
	rm a.out poly

clean:
	rm -rf src cpplot_*
	rm -f execname
