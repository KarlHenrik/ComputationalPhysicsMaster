#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "jacobiSolve.h"

using namespace  std;
using namespace  arma;

int main(int argc, char const *argv[]) {
    // Initializing matrices
    int n = 200;
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

    // Armadillo solution
    vec arma_eigval = zeros<vec>(n);
    mat arma_eigvecs = zeros<mat>(n, n);
    eig_sym(arma_eigval, arma_eigvecs, A);
    arma_eigval.save("Results/arma_eigval", arma_ascii);
    arma_eigvecs.save("Results/arma_eigvecs", arma_ascii);
    cout << "Armadillo solver: " << arma_eigval(0) << endl;

    // Diagonalizing with Jacobi rotations
    jacobiSolve(A, R, n);
    vec jacobi_eigvals_unsorted = diagvec(A);
    uvec eigvals_order = sort_index(jacobi_eigvals_unsorted);
    vec jacobi_eigvals = zeros<vec>(n);
    mat jacobi_eigvecs = zeros<mat>(n, n);
    for (int i = 0; i < n; i++) {
        jacobi_eigvals(i) = jacobi_eigvals_unsorted(eigvals_order(i));
        jacobi_eigvecs.col(i) = R.col(eigvals_order(i));
    }
    jacobi_eigvals.save("Results/jacobi_eigval", arma_ascii);
    jacobi_eigvecs.save("Results/jacobi_eigvecs", arma_ascii);
    cout << "Jacobi: " << jacobi_eigvals(0) << endl;

    // Analytical solution
    vec analytical_eigvals = zeros<vec>(n);
    mat analytical_eigvecs = zeros<mat>(n, n);
    double factor = M_PI / (n + 1);
    for (int col = 0; col < n; col++) {
        analytical_eigvals(col) = d + 2 * a * cos((col + 1) * factor);
        for (int row = 0; row < n; row++) {
            analytical_eigvecs(row, col) = sin((row + 1) * (col + 1) * factor);
        }
    }
    analytical_eigvals.save("Results/analytical_eigval", arma_ascii);
    analytical_eigvecs.save("Results/analytical_eigvecs", arma_ascii);
    cout << "Analytical: "  << d + 2 * a * cos(1 * M_PI / (n + 1)) << endl;
    return 0;
}
