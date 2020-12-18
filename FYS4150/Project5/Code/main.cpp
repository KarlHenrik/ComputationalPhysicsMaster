#include <vector>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "system.h"
#include "solver.h"

using namespace std;

int main(int argc, char const *argv[]) {
    /*
    {   // 1D simple model
        vector <int> steps_to_test {11, 101};
        for (int n : steps_to_test) {
            vec initial = zeros<vec>(n);
            initial(n - 1) = 1;
            double dx = 1.0 / (n - 1);
            double dt = dx * dx / 2;

            int T_001 = 1/dt / 100; // the number of iterations up to t = 0.01
            vector <int> intervals {T_001, T_001 * 9};
            vector <string> schemes {"exp", "imp", "crank"};

            for (string scheme : schemes) {
                System system(initial, dx, "Output/rod_" + scheme + "_" + to_string(n) + ".txt");
                solve(system, scheme, intervals, dt);
            }
        }
    }
    */

    {   // 2D simple model
        vector <int> steps_to_test {11, 101};
        for (int n : steps_to_test) {
            double dx = 1.0 / (n - 1);
            mat initial = zeros<mat>(n, n);
            for (int x = 0; x < n; x++) {
                for (int y = 0; y < n; y++) {
                    initial(x, y) = sin(M_PI * x * dx) * sin(M_PI * y * dx);
                }
            }
            double dt = dx * dx / 4.0;
            int T_001 = 1/dt / 100; // the number of iterations up to t = 0.01
            vector <int> intervals {T_001, T_001 * 4};
            vector <string> schemes {"exp", "imp"};

            for (string scheme : schemes) {
                System system(initial, dx, "Output/plate_" + scheme + "_" + to_string(n) + ".txt");
                solve(system, scheme, intervals, dt);
            }
        }
    }

    /*
    {   // Lithosphere model
        int n = 1000;
        vec initial = linspace<vec>(0, 1, n);
        double dx = 1.0 / n;
        double dt = 1.0 / 1000.0;

        vector <int> intervals {10, 490, 500, 9000};
        vector <string> sources {"base", "enrichedConst", "enriched"};
        double D = 1.5653;
        for (string source : sources) {
            System system(initial, dx, "Output/lit_" + source + ".txt", 0, D, source);
            solve(system, "crank", intervals, dt);
        }
    }
    */


    cout << "Finished sucessfully" << endl;
    return 0;
}
