LIB_DIR = lib

default: spinOps

spinOps: setup.py spinOps.pyx $(LIB_DIR)/libspinOps.a
	python3 setup.py build_ext --inplace && rm -f spinOps.c && rm -Rf build

$(LIB_DIR)/libspinOps.a:
	make -C $(LIB_DIR) libspinOps.a

clean:
	rm *.so
