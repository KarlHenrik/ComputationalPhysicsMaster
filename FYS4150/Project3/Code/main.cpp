#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarSystem.h"
#include "euler.h"
#include "vVerlet.h"
using namespace std;

int main(int argc, char const *argv[]) {
    int n = 500000;

    SolarSystem solarSystem("Output/ssys.txt");
    //solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1988500, "Sun");
    //solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219, "Earth");
    solarSystem.addFromFile("planetData.txt");

    //Euler solver(5.0 / n);
    vVerlet solver(250.0 / n);

    for (int i = 0; i < n; i++) {
        solarSystem.writePosVel();
        solver.step(solarSystem);
    }

    return 0;
}
