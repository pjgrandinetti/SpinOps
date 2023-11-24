CC = gcc

objects =  spinOps.o spatialOps.o

libspinOps.a : $(objects)
	ar rcs libspinOps.a $(objects)

spatialOps.o : spatialOps.h
	CC -c -O3 spatialOps.c

spinOps.o : spinOps.h
	CC -c -O3 spinOps.c

clean :
	rm $(objects) libspinOps.a 
