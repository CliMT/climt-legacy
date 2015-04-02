
CYTHON = _gfs_dynamics.pyx

CYTHONOBJS= _gfs_dynamics.c\
            _gfs_dynamics.so

cython_wrapper: $(CYTHON) 
	cd src; make python_lib
	LDFLAGS="-L"$(PWD)\/src python setup.py build_ext --inplace
	rm -R build/

all: $(cython_wrapper)

clean: 
	cd src; make clean
	rm $(CYTHONOBJS)

