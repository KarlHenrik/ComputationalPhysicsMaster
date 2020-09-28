#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;

void jacobiSolve(mat &, mat &, int);
void offdiag(mat &, int *, int *, int);
void Jacobi_rotate (mat &, mat &, int, int, int);
