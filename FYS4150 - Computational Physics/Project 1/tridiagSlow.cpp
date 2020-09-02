/* FYS4150 Project 1 - 1.b)
This program takes one additional input "k" which determines the numer of points n=10^k. It writes the results to a file.
Solving the one-dimensional Poisson equation with Dirichlet boundary.
This inefficient solution to the problem does not make use of the fact the diagonals
have the same predetermined values.
https://compphysics.github.io/ComputationalPhysics/doc/Projects/2020/Project1/html/Project1.html
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

    // Writing to file
    string ofilename = "Results/sl";
    ofilename += to_string(exponent);
    ofilename += ".txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    // x     u_solution       u_analytical      log10(relErr)
    ofile << setw(15) << "x";
    ofile << setw(15) << "u_solution";
    ofile << setw(15) << "u_analytical";
    ofile << setw(15) << "log10(relErr)" << endl;
    for (int i = 1; i < n; i++) { //Not including endpoints, as we don't want to divide by 0
        double xi = x[i];
        double relErr = fabs((exact(xi) - sol[i]) / exact(xi));

        ofile << setw(15) << setprecision(8) << xi;
        ofile << setw(15) << setprecision(8) << sol[i];
        ofile << setw(15) << setprecision(8) << exact(xi);
        ofile << setw(15) << setprecision(8) << log10(relErr) << endl;
    }
    ofile.close();

    delete [] d; delete [] er; delete [] el; delete [] b; delete [] x; delete [] sol;
    cout << "Successfully wrote results to: " << ofilename << endl;
    return 0;
}
