#ifndef Solver_H
#define Solver_H

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "system.h"

using namespace std;
using namespace arma;

void solve(System &, string, vector<int>, double);
void schemeChooser(System &, string, double, int);
void exp(System &, double, int);
void imp(System &, double, int);
void crankN(System &, double, int);
void exp2D(System &, double, int);
void imp2D(System &, double, int);


#endif // Solver_H
