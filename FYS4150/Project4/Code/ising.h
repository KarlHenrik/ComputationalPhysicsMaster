#ifndef Ising_H
#define Ising_H

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>

using namespace std;
using namespace arma;

class Ising {
    public:
        mat lattice;
        vec deltaEProb;
        int L, L2, A;
        double T, E, M, Arate;
        vec EVs = zeros<vec>(5);

        Ising(int, double, bool);
        void flipCycle(mt19937_64 &, uniform_real_distribution<double> &);
        void calcEVs(int);
        void writeEVs(ofstream &);
        void calcWriteChange(int, string);
        inline int periodic(int, int);

};

#endif // Ising_H
