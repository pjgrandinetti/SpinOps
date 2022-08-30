CC = gcc

objects =  spinOp.o

libspinOp.a : $(objects)
	ar rcs libspinOp.a $(objects)

spinOp.o : spinOp.h
	CC -c -O3 spinOp.c

clean :
	rm $(objects) spinOp.a 
