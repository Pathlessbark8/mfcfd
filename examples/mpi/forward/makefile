ifndef MF_DIR
	$(error directory not set)
endif

SRCDIR=$(MF_DIR)/src_forward

all:
	@echo "make explicit or implicit"

explicit: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp explicit/* . && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built explicit solver"

implicit: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && cp implicit/* . && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built explicit solver"

clean:
	rm -rf src
	rm -f execname
