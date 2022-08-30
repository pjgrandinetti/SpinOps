LIB_DIR = lib

default: pyexamples

pyexamples: setup.py pySpinOp.pyx $(LIB_DIR)/libspinOp.a
	python3 setup.py build_ext --inplace && rm -f pySpinOp.c && rm -Rf build

$(LIB_DIR)/libspinOp.a:
	make -C $(LIB_DIR) libspinOp.a

clean:
	rm *.so
