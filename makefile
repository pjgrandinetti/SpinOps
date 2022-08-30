LIB_DIR = lib

default: spinOpy

spinOpy: setup.py spinOpy.pyx $(LIB_DIR)/libspinOp.a
	python3 setup.py build_ext --inplace && rm -f spinOpy.c && rm -Rf build

$(LIB_DIR)/libspinOp.a:
	make -C $(LIB_DIR) libspinOp.a

clean:
	rm *.so
