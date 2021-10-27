CFLAGS = -O3
CC = g++
executables = main testing
libs = -larmadillo -llapack -lblas

all: main testing

# main
main: main.o system.o solver.o
	$(CC) $(CFLAGS) -fopenmp -o main main.o system.o solver.o

main.o: main.cpp
	$(CC) $(CFLAGS) -fopenmp -c main.cpp $(libs)

system.o: system.cpp
	$(CC) $(CFLAGS) -c system.cpp $(libs)

solver.o: solver.cpp
	$(CC) $(CFLAGS) -fopenmp -c solver.cpp $(libs)

# Testing
tests-main.o: UnitTests/tests-main.cpp
	$(CC) $(CFLAGS) -c UnitTests/tests-main.cpp

testing: testing.o tests-main.o system.o solver.o
	$(CC) $(CFLAGS) -fopenmp -o testing testing.o tests-main.o system.o solver.o $(libs)

testing.o: UnitTests/testing.cpp
	$(CC) $(CFLAGS) -c UnitTests/testing.cpp

clean:
	rm *.o $(executables)
