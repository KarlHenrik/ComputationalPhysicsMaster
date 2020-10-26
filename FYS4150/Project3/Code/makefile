CFLAGS = -O3
CC = g++
executables = main earthSun

all: main earthSun testing

# main
main: main.o vec3.o solarSystem.o celestialBody.o vVerlet.o
	$(CC) $(CFLAGS) -o main main.o vec3.o solarSystem.o celestialBody.o vVerlet.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

celestialBody.o: celestialBody.cpp
	$(CC) $(CFLAGS) -c celestialBody.cpp

solarSystem.o: solarSystem.cpp
	$(CC) $(CFLAGS) -c solarSystem.cpp

vVerlet.o: vVerlet.cpp
	$(CC) $(CFLAGS) -c vVerlet.cpp

euler.o: euler.cpp
	$(CC) $(CFLAGS) -c euler.cpp

# Utility
vec3.o: vec3.cpp
	$(CC) $(CFLAGS) -c vec3.cpp

# NonOO
earthSun: earthSun.o vec3.o
	$(CC) $(CFLAGS) -o earthSun earthSun.o vec3.o

earthSun.o: earthSun.cpp
	$(CC) $(CFLAGS) -c earthSun.cpp

# Testing
tests-main.o: UnitTests/tests-main.cpp
	$(CC) $(CFLAGS) -c UnitTests/tests-main.cpp

testing: testing.o tests-main.o vec3.o solarSystem.o celestialBody.o vVerlet.o euler.o
	$(CC) $(CFLAGS) -o testing testing.o tests-main.o vec3.o solarSystem.o celestialBody.o vVerlet.o euler.o

testing.o: UnitTests/testing.cpp
	$(CC) $(CFLAGS) -c UnitTests/testing.cpp

clean:
	rm *.o $(executables)
