all: buckling bucklingTiming oneElectron twoElectron oneElectronRhoN

buckling: buckling.o jacobiSolve.o
	g++ -O3 -o buckling buckling.o jacobiSolve.o -larmadillo -llapack -lblas

buckling.o: buckling.cpp
	g++ -O3 -c buckling.cpp

bucklingTiming: bucklingTiming.o jacobiSolve.o
	g++ -O3 -o bucklingTiming bucklingTiming.o jacobiSolve.o -larmadillo -llapack -lblas

bucklingTiming.o: bucklingTiming.cpp
	g++ -O3 -c bucklingTiming.cpp

oneElectron: oneElectron.o
	g++ -O3 -o oneElectron oneElectron.o jacobiSolve.o -larmadillo -llapack -lblas

oneElectron.o: oneElectron.cpp
	g++ -O3 -c oneElectron.cpp

twoElectron: twoElectron.o
	g++ -O3 -o twoElectron twoElectron.o jacobiSolve.o -larmadillo -llapack -lblas

twoElectron.o: twoElectron.cpp
	g++ -O3 -c twoElectron.cpp

oneElectronRhoN: oneElectronRhoN.o
	g++ -O3 -o oneElectronRhoN oneElectronRhoN.o jacobiSolve.o -larmadillo -llapack -lblas

oneElectronRhoN.o: oneElectronRhoN.cpp
	g++ -O3 -c oneElectronRhoN.cpp

jacobiSolve.o: jacobiSolve.cpp
	g++ -O3 -c jacobiSolve.cpp

clean:
	rm *.o buckling bucklingTiming oneElectron twoElectron oneElectronRhoN
