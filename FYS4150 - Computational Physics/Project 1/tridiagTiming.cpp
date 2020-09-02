/* FYS4150 Project 1 - 1.c) + 1.d))
This program takes one additional input "k" which determines the numer of points n=10^k. It writes the results to a file.
Solving the one-dimensional Poisson equation with Dirichlet boundary.
This program uses the predetermined entries in the tridiagonal matrix to speed up the algorithm.
https://compphysics.github.io/ComputationalPhysics/doc/Projects/2020/Project1/html/Project1.html
c++ -o proj1b.exe proj1b.cpp && proj1b 1bout 1 && proj1b 1bout 2 && proj1b 1bout 3
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "lib.h"

using namespace std;

double tridiag(int exponent);
double tridiagSlow(int exponent);
double tridiagLU(int exponent);

double f(double x) {
    return 100 * exp(-10 * x);
}

double exact(double x) {
    return 1 - (1 - exp(-10)) * x - exp(-10 * x);
}

int main(int argc, char const *argv[]) {
    // Writing to file
    string ofilename = "Results/timing.txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    // x     u solution       u analytical      log10(relErr)
    ofile << setw(15) << "n";
    ofile << setw(15) << "tridiag";
    ofile << setw(15) << "tridiagSlow";
    ofile << setw(15) << "tridiagLU" << endl;
    for (int i = 1; i < 11; i++) { //Not including endpoints, as we don't want to divide by 0
        ofile << setw(15) << setprecision(8) << i;
        ofile << setw(15) << setprecision(8) << tridiag(i);
        ofile << setw(15) << setprecision(8) << tridiagSlow(i);
        if (i < 5) {
            ofile << setw(15) << setprecision(8) << tridiagLU(i);
        }
        ofile << endl;
    }
    ofile.close();
    cout << "Successfully wrote results to: " << ofilename << endl;
    return 0;
}

double tridiag(int exponent) {
    // constants
    int n = pow(10.0, exponent); //points
    double h = 1.0 / n; //steplength
    double hh = h * h;
    // arrays
    double *d = new double[n+1]; //main diagonal
    double *b = new double[n+1]; //right side
    double *x = new double[n+1]; //x-values
    double *sol = new double[n+1]; //answer
    // Initial values
    for (int i = 1; i < n + 1; i++) { // the endpoints are never used (index 0 and n+1)
        x[i] = i * h; // x is not used in the implementation of the algorithm
        b[i] = hh * f(i * h);
        d[i] = (i + 1.0) / i; // d is pre-calculated according to a formula
    }
    clock_t start, finish;
    start = clock();
    // Forward substitution. Updating all diagonal elements and right side elements, from the second row
    for (int i = 2; i < n; i++) {
        b[i] = b[i] + b[i - 1] / d[i - 1];
    }
    // Backward substitution. With only two diagonals, we can find the coefficients moving backwards
    sol[n-1] = b[n-1] / d[n-1];
    for (int i = n - 2; i > 0; i--) {
        sol[i] = (b[i] + sol[i + 1]) / d[i];
    }
    finish = clock();
    double timeused = (finish - start) / ((double) CLOCKS_PER_SEC);
    delete [] d; delete [] b; delete [] x; delete [] sol;

    return timeused;
}

double tridiagSlow(int exponent) {
    // constants
    int n = pow(10, exponent); //points
    double h = 1.0 / n; //steplength
    double hh = h * h;
    // arrays
    double *d = new double[n+1]; //main diagonal
    double *er = new double[n+1]; //left diagonal
    double *el = new double[n+1]; //right diagonal
    double *b = new double[n+1]; //right side
    double *x = new double[n+1]; //x-values
    double *sol = new double[n+1]; //answer
    // Initial values
    sol[0] = 0;
    sol[n] = 0;
    for (int i = 1; i < n; i++) { // d and e have elements from 1 to n (indexes 1 to n-1)
        d[i] = 2;
        er[i] = -1;
        el[i] = -1;
    }
    for (int i = 0; i < n + 1; i++) {
        x[i] = i * h; // x is not used in the implementation of the algorithm
        b[i] = hh * f(i * h);
    }
    clock_t start, finish;
    start = clock();
    // Forward substitution. Updating all diagonal elements and right side elements, from the second row
    double factor;
    for (int i = 2; i < n; i++) {
        factor = el[i - 1] / d[i - 1]; //the factor used for the elimination at each step
        d[i] = d[i] - er[i - 1] * factor;
        b[i] = b[i] - b[i - 1] * factor;
    }
    // Backward substitution. With only two diagonals, we can find the coefficients moving backwards
    sol[n-1] = b[n-1] / d[n-1];
    for (int i = n - 2; i > 0; i--) {
        sol[i] = (b[i] - er[i] * sol[i + 1]) / d[i];
    }
    finish = clock();
    double timeused = (finish - start) / ((double) CLOCKS_PER_SEC);

    delete [] d; delete [] er; delete [] el; delete [] b; delete [] x; delete [] sol;
    return timeused;
}

double tridiagLU(int exponent) {
    // constants
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
    clock_t start, finish;
    start = clock();
    // Solving
    int *indx = new int[n]; double d; // return variables we won't use
    ludcmp(A, n, indx, &d); // LU decompose  a[][]
    lubksb(A, n, indx, b); // Solving equations
    finish = clock();
    double timeused = (finish - start) / ((double) CLOCKS_PER_SEC);
    // Deleting arrays
    for (int i = 0; i < n; i++) {
        delete [] A[i];
    }
    delete [] A; delete [] b; delete [] x; delete [] indx;

    return timeused;
}
