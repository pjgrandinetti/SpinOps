CC = gcc

objects =  spinOps.o

libspinOps.a : $(objects)
	ar rcs libspinOps.a $(objects)

spinOps.o : spinOps.h
	CC -c -O3 spinOps.c

clean :
	rm $(objects) libspinOps.a 
