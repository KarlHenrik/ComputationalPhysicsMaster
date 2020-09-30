#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

#include "jacobiSolve.h"

TEST_CASE( "Finding largest element of a symmetric matrix", "[offdiag]" ) {
    int p; int q;
    // Initializing matrix
    int n = 3;
    mat A = zeros<mat>(n,n);
    A(0,0) = 100;
    A(1,0) = 5;
    A(0,1) = 5;
    A(1,2) = -7; // Negative but still "largest"
    A(2,1) = -7;

    offdiag(A, &p, &q, n);
    REQUIRE( p == 1 );
    REQUIRE( q == 2 );
    // New largest
    A(1,0) = 8;
    A(0,1) = 8;
    offdiag(A, &p, &q, n);
    REQUIRE( p == 0 );
    REQUIRE( q == 1 );
}

TEST_CASE( "Single rotations sets elements to zero", "[Jacobi_rotate]" ) {
    int p; int q;
    // Initializing matrix
    int n = 3;
    mat A = zeros<mat>(n,n);
    mat R = eye(n, n);
    A(0,0) = 100;
    A(1,0) = 5;
    A(0,1) = 5;
    A(1,2) = -7; // Negative but still "largest"
    A(2,1) = -7;

    offdiag(A, &p, &q, n);
    Jacobi_rotate(A, R, p, q, n);
    REQUIRE( A(p, q) == 0 );
    REQUIRE( A(q, p) == 0 );
}

TEST_CASE( "First three buckling beam eigvals using Jacobi", "[jacobiSolve]" ) {
    int n = 100;
    mat R = eye(n, n); // The matrix R which will hold all of the eigenvectors
    mat A = zeros<mat>(n,n); // The tridiagonal matrix A which is to be diagonalized
    double h = 1.0 / (n + 1);
    double h2 = h * h;
    double d = 2 / h2;
    double a = - 1 / h2;
    A(0, 0) = d;
    for (int row = 1; row < n; row++) {
        A(row, row) = d;
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
    REQUIRE( jacobi_eigvals(0) == Approx(9.86880867886066) );
    REQUIRE( jacobi_eigvals(1) == Approx(39.4656872804080) );
    REQUIRE( jacobi_eigvals(2) == Approx(88.7620027360826) );

}

TEST_CASE( "Orthogonality of eigenvectors for one electron", "[jacobiSolve]" ) {
    // Constants
    int n = 300;
    double rho_N = 6;
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
    mat jacobi_eigvecs = zeros<mat>(n, n);
    for (int i = 0; i < n; i++) {
        jacobi_eigvals(i) = jacobi_eigvals_unsorted(eigvals_order(i));
        jacobi_eigvecs.col(i) = R.col(eigvals_order(i));
    }
    REQUIRE( jacobi_eigvecs(0) * jacobi_eigvecs(1) == Approx(0).epsilon(0.001) );
    REQUIRE( jacobi_eigvecs(23) * jacobi_eigvecs(55) == Approx(0).epsilon(0.001) );
    REQUIRE( jacobi_eigvecs(201) * jacobi_eigvecs(14) == Approx(0).epsilon(0.001) );
}
