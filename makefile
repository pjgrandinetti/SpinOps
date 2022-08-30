CC = gcc

objects =  spinOp.o

libexamples.a : $(objects)
	ar rcs spinOp.a $(objects)

spinOp.o : spinOp.h
	CC -c -O3 spinOp.c

clean :
	rm $(objects) spinOp.a 
