#include <vector>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include "ising.h"
#include "tempTester.h"
#include "omp.h"
#include "time.h"

using namespace std;

int main(int argc, char const *argv[]) {


    // Testing 2x2 model
    TempTester tester2x2(2, 1, 4 + 0.0000001, 0.01, true); //(int L, double T0, double TN, double dT, bool random)
    tester2x2.calcParallell(100000, 12);
    tester2x2.write("Output/2x2.txt");


    // Timing
    for (int L = 20; L <= 100; L += 20) {
        TempTester testerTimeParalell(L, 2.2, 2.4 + 0.0000001, 0.05, true); //(int L, double T0, double TN, double dT, bool random)
        cout << "Timing parallell with L = " << L << " ";
        testerTimeParalell.calcParallell(100000, 12);
        TempTester testerTimeSingle(L, 2.2, 2.4 + 0.0000001, 0.05, true); //(int L, double T0, double TN, double dT, bool random)
        cout << "Timing non-parallell with L = " << L << " ";
        testerTimeSingle.calcParallell(100000, 1);
    }


    // Testing 20x20 to 100x100 models - THIS TAKES A LONG TIME!!
    for (int L = 20; L <= 100; L += 20) {
        TempTester testerBig(L, 2.2, 2.4 + 0.0000001, 0.001, true); //(int L, double T0, double TN, double dT, bool random)
        testerBig.calcParallell(1000000, 12);
        testerBig.write("Output/" + to_string(L) + "e7rnd.txt");
    }


    // Finding time to equilibrium for 20x20
    Ising mdl1(20, 1, false); //(int L_, double T_, bool random)
    mdl1.calcWriteChange(100000, "Output/tte/T1.txt");
    Ising mdl2(20, 1, true);
    mdl2.calcWriteChange(100000, "Output/tte/T1rand.txt");
    Ising mdl3(20, 2.4, false);
    mdl3.calcWriteChange(100000, "Output/tte/T24.txt");
    Ising mdl4(20, 2.4, true);
    mdl4.calcWriteChange(100000, "Output/tte/T24rand.txt");

    cout << "Finished sucessfully" << endl;
    return 0;
}
