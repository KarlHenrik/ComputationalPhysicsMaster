# Project 5

## Code Folder

We have written code for solving the diffusion equation in one or two dimensions. In the main program, we define the initial and boundrary conditions of our models, and then create a `system.cpp` object to hold the values and boundary conditions. `system.cpp` also contains the different source functions for the lithosphere model, and functionality for writing to file which is called from `solver.cpp`.

We give the `system` object to the `solver.solve()` function together with the scheme we want to use, the time step and how many iterations to run the solver between each write to file. `solver.solve()` then takes care of the rest, leading to a hopefully compact and user friendly making and running of models in `main`.

`testing.cpp` contains unit tests for `system.cpp` and `solver.cpp`.

The jupyter notebooks read the results from the c++ calculations from files and plots them and compares them to analytical solutions.

You can compile the code by running the "make" command using the makefile. The `main` and `testing` programs are run with no input parameters.

## Report Folder

Contains our report pdf.

### About implementation choices
`solver.cpp` has similar code repeated a couple of times (the forward and backward substitution), but I preferred it this way as there is no `solver` object to hold all the different variables, and the "reused code" is so short that I felt it was more readable as is.
