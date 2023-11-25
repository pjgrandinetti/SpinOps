CC = gcc
CFLAGS = -I/opt/homebrew/Cellar/cunit/2.1-3/include
LDFLAGS = -L/opt/homebrew/lib -lcunit
ARFLAGS = rcs

# List of objects
OBJECTS =  spinOps.o spatialOps.o

# Name of the library
LIBRARY = libspinOps.a

# Default target
all: $(LIBRARY)

# Rule to create the library
$(LIBRARY) : $(OBJECTS)
	ar $(ARFLAGS) $(LIBRARY) $(OBJECTS)

# Rule to compile source files
%.o : %.c %.h
	$(CC) $(CFLAGS) -c -O3 $< -o $@

# Rule to clean the objects and the library
clean :
	rm -f $(OBJECTS) $(LIBRARY)

# Rule to compile and run the tests
test: test.c $(OBJECTS)
	$(CC) -o test test.c $(OBJECTS) $(CFLAGS) $(LDFLAGS)
