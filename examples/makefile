ifndef FF_DIR
	$(error directory not set)
endif

SRCDIR=$(FF_DIR)/src


all: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -j -f makefile execname 
	mv src/execname .
	@echo "Built solver"

clean:
	rm -rf src
	rm -f execname
