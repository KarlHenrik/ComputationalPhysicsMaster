#ifndef TempTester_H
#define TempTester_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include "ising.h"
#include "omp.h"

using namespace std;
using namespace arma;

class TempTester {
    public:
        vector<Ising> models;

        TempTester(int, double, double, double, bool);
        void calc(int);
        void calcParallell(int, int);
        void write(string);
};

#endif // TempTester_H
