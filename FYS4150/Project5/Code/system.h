#ifndef System_H
#define System_H

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace arma;

class System {
    public:
        mat U;
        vec u;
        double dx, t, D;
        int t_i, dims;
        ofstream ofile;
        string sourceName;
        double Qfac = 6.9785;// * pow(10, 7);

        System(vec &, double, string, double = 0, double = 1, string = "none");
        System(mat, double, string, double = 0, double = 1);
        void writeState();
        void close();

        double source(double, double);
        double heatProd(double);
        double heatProdEnrichedConst(double);
        double heatProdEnriched(double, double);
};

#endif // System_H
