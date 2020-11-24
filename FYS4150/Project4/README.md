# Project 4

The report can be found in the Report folder.

## Code

### C++ files
You compile the code by running the make command, in the folder with the makefile. The main code is run with no input parameters. The testing code is run with no input parameters.

main.cpp creates the Ising models of interest to the report and writes the results to files in the Output folder.

ising.cpp is a class that lets you create objects that can be thought of as instance of an Ising model grid, with a specific temperature. Each grid size and temperature requires its own ising object. The constructor decides the temperature, size and if the spins should be initialized randomly. By calling calcEVs(N), you perform N sweeps of the lattice, and update all expectation values for each sweep. writeEVs() then lets you write the results to a file. calcWriteChange() writes the values of the system to file for each sweep.

tempTester.cpp lets you create objects that hold many Ising objects. The constructor takes in the grid size and temperature RANGE to make models for. calcParallell() then calls calcEVs() on each object, assigning an equal number of Ising objects to each thread. calc() calls calcEVs for all objects without parallelization. writeEVs() then writes the results from all the models to a file in order of temperature.

The use cases in main.cpp should be self explanatory.

### Jupyter notebook
P4Notebook.ipynb plots the results from the c++ codes, and also does a little handling of the data, like finding the critical temperatures.
