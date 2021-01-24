CC=clang++
CFALGS=-Wall -std=c++17
LDFLAGS= -lgmpxx -lgmp
HEADERS = polynomial.h
OBJECTS = main.o polynomial.o

defualt: polyroot

%.o : %.cpp $(HEADERS)
	$(CC) -c $(CFLAGS) $< 

polyroot: $(OBJECTS)
	$(CC) $(CFALGS) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm -f *.o
	rm -f polyroot
