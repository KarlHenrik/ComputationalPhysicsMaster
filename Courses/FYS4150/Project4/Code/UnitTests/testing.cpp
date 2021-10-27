#include "catch.hpp"

#include "../tempTester.h"
#include "../ising.h"

TEST_CASE("Non-random lattice all ones", "[Ising(int, double, bool)]") {
    Ising mdl1(20, 1, false); //(int L_, double T_, bool random)
    int mistmatch = 0;
    for (int x = 0; x < 20; x++) {
        for (int y = 0; y < 20; y++) {
            if (mdl1.lattice(x, y) != 1) {
                mistmatch++;
            }
        }
    }
    REQUIRE( mistmatch == 0 );
}

TEST_CASE("Random lattice not all ones (could randomly fail)", "[Ising(int, double, bool)]") {
    Ising mdl1(100, 1, true); //(int L_, double T_, bool random)
    int mistmatch = 0;
    for (int x = 0; x < 20; x++) {
        for (int y = 0; y < 20; y++) {
            if (mdl1.lattice(x, y) != 1) {
                mistmatch++;
            }
        }
    }
    REQUIRE( mistmatch != 0 );
}

TEST_CASE("T = 1 leads to low energy (could randomly fail)", "[calcEVs(int)]") {
    Ising mdl1(20, 1, true); //(int L_, double T_, bool random)
    mdl1.calcEVs(3000);
    REQUIRE( mdl1.EVs(0) < -1.9 );
}

TEST_CASE("T = 3 leads to higher energy (could randomly fail)", "[calcEVs(int)]") {
    Ising mdl1(20, 3, true); //(int L_, double T_, bool random)
    mdl1.calcEVs(3000);
    REQUIRE( mdl1.EVs(0) > -1.5 );
}

TEST_CASE("T = 3 leads to low magnetization (could randomly fail)", "[calcEVs(int)]") {
    Ising mdl1(20, 3, true); //(int L_, double T_, bool random)
    mdl1.calcEVs(3000);
    REQUIRE( mdl1.EVs(4) < 0.3 );
}

TEST_CASE("Parallelized code leads to similar results as non-parallelized, and single thread (could randomly fail)", "[calcEVs(int)]") {
    TempTester testerPara(20, 2.2, 2.4 + 0.0000001, 0.05, true); //(int L, double T0, double TN, double dT, bool random)
    testerPara.calcParallell(1000, 12);
    TempTester testerParaOne(20, 2.2, 2.4 + 0.0000001, 0.05, true); //(int L, double T0, double TN, double dT, bool random)
    testerParaOne.calcParallell(1000, 1);
    TempTester tester(20, 2.2, 2.4 + 0.0000001, 0.05, true); //(int L, double T0, double TN, double dT, bool random)
    tester.calc(1000);
    REQUIRE( testerPara.models[0].EVs(0) == Approx(tester.models[0].EVs(0)).epsilon(0.1) );
    REQUIRE( testerPara.models[0].EVs(0) == Approx(testerParaOne.models[0].EVs(0)).epsilon(0.1) );
    REQUIRE( testerParaOne.models[0].EVs(0) == Approx(tester.models[0].EVs(0)).epsilon(0.1) );
}
