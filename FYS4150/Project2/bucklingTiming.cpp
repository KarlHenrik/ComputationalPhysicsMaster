#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "time.h"

#include "jacobiSolve.h"

using namespace  std;
using namespace  arma;

double armaTime(mat, int);
double jacobiTime(mat, mat, int);
double analyticalTime(int, double, double);

int main(int argc, char const *argv[]) {
    // Writing to file
    string ofilename = "Results/bucklingTiming.txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    // x     u solution       u analytical      log10(relErr)
    ofile << setw(15) << "n";
    ofile << setw(15) << "Armadillo";
    ofile << setw(15) << "Analytical";
    ofile << setw(15) << "Jacobi" << endl;
    for (int n = 100; n < 1000; n = n + 100) {
        mat R = eye(n, n); // The matrix R which will hold all of the eigenvectors
        mat A = zeros<mat>(n,n); // The tridiagonal matrix A which is to be diagonalized
        double h2 = 1.0 / (n * n);
        double d = 2 / h2;
        double a = - 1 / h2;
        A(0, 0) = d;
        for (int row = 1; row < n; row++) {
            A(row, row) = d;
            A(row - 1, row) = a;
            A(row, row - 1) = a;
        }

        ofile << setw(15) << setprecision(8) << n;
        ofile << setw(15) << setprecision(8) << armaTime(A, n);
        ofile << setw(15) << setprecision(8) << analyticalTime(n, d, a);
        if (n <= 400) {
            ofile << setw(15) << setprecision(8) << jacobiTime(A, R, n);
        }
        ofile << endl;
    }
    ofile.close();
    cout << "Successfully wrote results to: " << ofilename << endl;
    return 0;
}

double armaTime(mat A, int n) {
    vec arma_eigval = zeros<vec>(n);
    mat arma_eigvecs = zeros<mat>(n, n);
    clock_t start, finish;
    start = clock();
    eig_sym(arma_eigval, arma_eigvecs, A);
    finish = clock();
    double timeused = (finish - start) / ((double) CLOCKS_PER_SEC);
    return timeused;
}

double jacobiTime(mat A, mat R, int n) {
    clock_t start, finish;
    start = clock();
    jacobiSolve(A, R, n);
    finish = clock();
    double timeused = (finish - start) / ((double) CLOCKS_PER_SEC);
    return timeused;
}

double analyticalTime(int n, double d, double a) {
    vec analytical_eigvals = zeros<vec>(n);
    mat analytical_eigvecs = zeros<mat>(n, n);
    double factor = M_PI / (n + 1);
    for (int col = 0; col < n; col++) {
        analytical_eigvals(col) = d + 2 * a * cos((col + 1) * factor);
        for (int row = 0; row < n; row++) {
            analytical_eigvecs(row, col) = sin((row + 1) * (col + 1) * factor);
        }
    }
}
