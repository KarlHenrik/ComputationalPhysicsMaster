/* FYS4150 Project 1 - 1.c) + 1.d))
This program takes one additional input "k" which determines the numer of points n=10^k. It writes the results to a file.
Solving the one-dimensional Poisson equation with Dirichlet boundary.
This program uses functions from the lib.cpp program to solve the equations using an LU factorization.
https://compphysics.github.io/ComputationalPhysics/doc/Projects/2020/Project1/html/Project1.html
c++ -o proj1b.exe proj1b.cpp && proj1b 1bout 1 && proj1b 1bout 2 && proj1b 1bout 3
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "lib.h"

using namespace std;

double f(double x) {
    return 100 * exp(-10 * x);
}

double exact(double x) {
    return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        cout << "Missing n-size (10^k exponent)" << endl;
        exit(1);
    }
    // constants
    int exponent = atoi(argv[1]);
    int n = pow(10.0, exponent); //points
    double h = 1.0 / n; //steplength
    double hh = h * h;
    n = n - 1; //shifting max index 1 down to ignore x=1
    // Initializing
    double **A = new double*[n]; //nxn matrix
    double *b = new double[n]; //right side
    double *x = new double[n]; //x-values
    // Initial values
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        b[i] = hh * f((i + 1) * h); //shifting index 1 up to ignore x=0
        x[i] = (i + 1) * h;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0;
        }
    }
    for (int i = 0; i < n-1; i++) { // Filling in values of the tridiagonal A
        A[i][i] = 2;
        A[i][i+1] = -1;
        A[i+1][i] = -1;
    }
    A[n-1][n-1] = 2;
    // Solving
    int *indx = new int[n]; double d; // return variables we won't use
    ludcmp(A, n, indx, &d); // LU decompose  a[][]
    lubksb(A, n, indx, b); // Solving equations
    // Writing to file
    string ofilename = "Results/lu";
    ofilename += to_string(exponent);
    ofilename += ".txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    // x     u solution       u analytical      log10(relErr)
    ofile << setw(15) << "x";
    ofile << setw(15) << "u_solution";
    ofile << setw(15) << "u_analytical";
    ofile << setw(15) << "log10(relErr)" << endl;
    for (int i = 0; i < n; i++) { //x = 0 and x = 1 are not included by default
        double xi = x[i];
        double relErr = fabs((exact(xi) - b[i]) / exact(xi));

        ofile << setw(15) << setprecision(8) << xi;
        ofile << setw(15) << setprecision(8) << b[i];
        ofile << setw(15) << setprecision(8) << exact(xi);
        ofile << setw(15) << setprecision(8) << log10(relErr) << endl;
    }
    ofile.close();
    // Deleting arrays
    for (int i = 0; i < n; i++) {
        delete [] A[i];
    }
    delete [] A; delete [] b; delete [] x; delete [] indx;
    cout << "Successfully wrote results to: " << ofilename << endl;
    return 0;
}
