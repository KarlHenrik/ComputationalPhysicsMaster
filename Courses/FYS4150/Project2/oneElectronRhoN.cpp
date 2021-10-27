#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "jacobiSolve.h"

using namespace  std;
using namespace  arma;

int main(int argc, char const *argv[]) {
    string ofilename = "Results/oneElectronRhoN.txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << "n/rho_N";
    ofile << setw(15) << "3";
    ofile << setw(15) << "6";
    ofile << setw(15) << "9";
    ofile << setw(15) << "12" << endl;
    for (int n = 50; n <= 250; n += 50) {
        ofile << setw(15) << setprecision(8) << n;
        for (double rho_N = 3; rho_N <= 12; rho_N += 3) {
            // Initializing matrices
            mat R = eye(n, n); // The matrix R which will hold all of the eigenvectors
            mat A = zeros<mat>(n,n); // The tridiagonal matrix A which is to be diagonalized
            double h = rho_N / (n + 1);
            double h2 = h * h;
            double d = 2 / h2;
            double a = -1 / h2;
            A(0, 0) = d + h2;
            for (int row = 1; row < n; row++) {
                double rho = (row + 1) * h;
                A(row, row) = d + rho * rho;
                A(row - 1, row) = a;
                A(row, row - 1) = a;
            }

            // Diagonalizing with Jacobi rotations
            jacobiSolve(A, R, n);
            vec jacobi_eigvals_unsorted = diagvec(A);
            uvec eigvals_order = sort_index(jacobi_eigvals_unsorted);
            vec jacobi_eigvals = zeros<vec>(n);
            for (int i = 0; i < n; i++) {
                jacobi_eigvals(i) = jacobi_eigvals_unsorted(eigvals_order(i));
            }
            ofile << setw(15) << setprecision(8) << fabs(3 - jacobi_eigvals(0));
            //cout << "n = " << n << " rho_N = " << rho_N << " lmda = " << fabs(3 - jacobi_eigvals(0)) << endl;
        }
        ofile << endl;
    }

    ofile.close();
    cout << "Successfully wrote results to: " << ofilename << endl;

    return 0;
}
