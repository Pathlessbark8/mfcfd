ifndef FF_DIR
	$(error directory not set)
endif

SRCDIR=$(FF_DIR)/src_mpi








all: 
	rm -rf src
	cp -r $(SRCDIR) src
	cd src && make -f makefile clean
	cp parameter.F90 src/
	cd src && make -f makefile execname 
	mv src/execname .

clean:
	rm -rf src
	rm -f execname
