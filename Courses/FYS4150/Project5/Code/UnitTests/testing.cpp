#include "catch.hpp"

#include "../system.h"
#include "../solver.h"

TEST_CASE("1D system converges to linear with all schemes", "[system, solve()]") {
    int n = 101;
    vec initial = zeros<vec>(n);
    initial(n - 1) = 1;
    double dx = 1.0 / (n - 1);
    double dt = dx * dx / 2;

    vector <int> intervals {20000}; // iterates to t = 1, where it should be pretty linear
    vector <string> schemes {"exp", "imp", "crank"};

    for (string scheme : schemes) {
        System system(initial, dx, "Output/test.txt");
        solve(system, scheme, intervals, dt);
        REQUIRE( system.u(25) == Approx(0.25).epsilon(0.001) );
        REQUIRE( system.u(50) == Approx(0.5).epsilon(0.001) );
        REQUIRE( system.u(75) == Approx(0.75).epsilon(0.001) );
    }
}

TEST_CASE("2D system converges to 0 schemes", "[system, solve()]") {
    int n = 11;
    double dx = 1.0 / (n - 1);
    mat initial = zeros<mat>(n, n);
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            initial(x, y) = sin(M_PI * x * dx) * sin(M_PI * y * dx);
        }
    }
    double dt = dx * dx / 4.0;
    int T = 1/dt; // the number of iterations up to t = 1
    vector <int> intervals {T};
    vector <string> schemes {"exp", "imp"};

    for (string scheme : schemes) {
        System system(initial, dx, "Output/test.txt");
        solve(system, scheme, intervals, dt);
        REQUIRE( system.U(2,2)+1 == Approx(1).epsilon(0.00001) );
        REQUIRE( system.U(4,6)+1 == Approx(1).epsilon(0.00001) );
        REQUIRE( system.U(5,5)+1 == Approx(1).epsilon(0.00001) );
    }
}

TEST_CASE("Source terms working as intended", "[system, source()]") {
    int n = 11;
    vec initial = zeros<vec>(n);
    initial(n - 1) = 1;
    double dx = 1.0 / (n - 1);
    double dt = dx * dx / 2;

    {
        System system(initial, dx, "Output/test.txt", 0, 1, "base");
        REQUIRE( system.source(0, 0) == Approx(1.4 * system.Qfac).epsilon(0.00001) );
        REQUIRE( system.source(1, 0) == Approx(0.05 * system.Qfac).epsilon(0.00001) );
    }
    {
        System system(initial, dx, "Output/test.txt", 0, 1, "enrichedConst");
        REQUIRE( system.source(0, 0) == Approx(1.4 * system.Qfac).epsilon(0.00001) );
        REQUIRE( system.source(1, 0) == Approx(0.55 * system.Qfac).epsilon(0.00001) );
    }
    {
        System system(initial, dx, "Output/test.txt", 0, 1, "enriched");
        REQUIRE( system.source(0.5, 300) == Approx(0.05 * system.Qfac).epsilon(0.00001) ); // enrichement should be gone
        REQUIRE( system.source(0.5, 0) != Approx(0.05 * system.Qfac).epsilon(0.00001) ); // enrichement should not be gone
    }
}
