CFLAGS = -O3
CC = g++
executables = main testing
libs = -larmadillo -llapack -lblas

all: main testing

# main
main: main.o ising.o tempTester.o #mpiexec -n ./main
	$(CC) $(CFLAGS) -fopenmp -o main main.o ising.o tempTester.o

main.o: main.cpp
	$(CC) $(CFLAGS) -fopenmp -c main.cpp $(libs)

ising.o: ising.cpp
	$(CC) $(CFLAGS) -c ising.cpp $(libs)

tempTester.o: tempTester.cpp
	$(CC) $(CFLAGS) -fopenmp -c tempTester.cpp

# Testing
tests-main.o: UnitTests/tests-main.cpp
	$(CC) $(CFLAGS) -c UnitTests/tests-main.cpp

testing: testing.o tests-main.o ising.o tempTester.o
	$(CC) $(CFLAGS) -fopenmp -o testing testing.o tests-main.o ising.o tempTester.o $(libs)

testing.o: UnitTests/testing.cpp
	$(CC) $(CFLAGS) -c UnitTests/testing.cpp

clean:
	rm *.o $(executables)
