#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>

#include "jacobiSolve.h"

using namespace  std;
using namespace  arma;

void twoElectronOmega(double);

int main(int argc, char const *argv[]) {
    // Constants
    int n = 200;
    double rho_N = 10;
    double omega = 1;
    // Initializing matrices
    mat R = eye(n, n); // The matrix R which will hold all of the eigenvectors
    mat A = zeros<mat>(n,n); // The tridiagonal matrix A which is to be diagonalized
    double h = rho_N / (n + 1);
    double h2 = h * h;
    double d = 2 / h2;
    double a = - 1 / h2;
    A(0, 0) = d + h2;
    for (int row = 1; row < n; row++) {
        double rho = (row + 1) * h;
        A(row, row) = d + omega * rho * rho + 1 / rho;
        A(row - 1, row) = a;
        A(row, row - 1) = a;
    }

    // Armadillo solution
    vec arma_eigval = zeros<vec>(n);
    mat arma_eigvecs = zeros<mat>(n, n);
    eig_sym(arma_eigval, arma_eigvecs, A);
    arma_eigval.save("Results/2Earma_eigval", arma_ascii);
    arma_eigvecs.save("Results/2Earma_eigvecs", arma_ascii);
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
    jacobi_eigvals.save("Results/2Ejacobi_eigval", arma_ascii);
    jacobi_eigvecs.save("Results/2Ejacobi_eigvecs", arma_ascii);
    cout << "Jacobi: " << jacobi_eigvals(0) << endl;

    twoElectronOmega(0.01);
    twoElectronOmega(0.5);
    twoElectronOmega(1);
    twoElectronOmega(5);

    return 0;
}

void twoElectronOmega(double omega) {
    // Constants
    int n = 200;
    double rho_N = 10;
    // Initializing matrices
    mat R = eye(n, n); // The matrix R which will hold all of the eigenvectors
    mat A = zeros<mat>(n,n); // The tridiagonal matrix A which is to be diagonalized
    double h = rho_N / (n + 1);
    double h2 = h * h;
    double d = 2 / h2;
    double a = - 1 / h2;
    A(0, 0) = d + h2;
    for (int row = 1; row < n; row++) {
        double rho = (row + 1) * h;
        A(row, row) = d + omega * rho * rho + 1 / rho;
        A(row - 1, row) = a;
        A(row, row - 1) = a;
    }
    // Diagonalizing with Jacobi rotations
    jacobiSolve(A, R, n);
    // Saving the ground state vector to file
    vec jacobi_eigvals_unsorted = diagvec(A);
    uvec eigvals_order = sort_index(jacobi_eigvals_unsorted);
    vec gs = R.col(eigvals_order(0));
    string filename = "Results/2EGS" + to_string(omega);
    gs.save(filename, arma_ascii);
}
