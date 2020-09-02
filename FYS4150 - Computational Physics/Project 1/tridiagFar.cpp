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
    // Forward substitution. Updating all diagonal elements and right side elements, from the second row
    for (int i = 2; i < n; i++) {
        b[i] = b[i] + b[i - 1] / d[i - 1];
    }
    // Backward substitution. With only two diagonals, we can find the coefficients moving backwards
    sol[n-1] = b[n-1] / d[n-1];
    for (int i = n - 2; i > 0; i--) {
        sol[i] = (b[i] + sol[i + 1]) / d[i];
    }
    // Writing to file
    double relErr = fabs((exact(x[n / 10]) - sol[n / 10]) / exact(x[n / 10]));
    std::cout << relErr << '\n';

    delete [] d; delete [] b; delete [] x; delete [] sol;
    cout << "Successfully wrote results to: " << endl;
    return 0;
}
