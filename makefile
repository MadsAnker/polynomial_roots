CC=clang++
HEADERS = polynomial.h

defualt: polyroot

polynomial.o: polynomial.cpp $(HEADERS)
	$(CC) -c polynomial.cpp -std=c++17 -o polynomial.o

main.o: main.cpp
	$(CC) -c main.cpp -std=c++17 -o main.o

polyroot: main.o polynomial.o
	$(CC) main.o polynomial.o -std=c++17 -lgmpxx -lgmp -o polyroot

clean:
	rm -f *.o
	rm -f polyroot
