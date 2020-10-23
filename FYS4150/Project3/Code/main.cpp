#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarSystem.h"
#include "euler.h"
#include "vVerlet.h"
using namespace std;

int main(int argc, char const *argv[]) {
    int n = 100000;
    double totalTime = 100;
    double dt = totalTime / n;
    double solarMass = 1988500;

    SolarSystem solarSystem("Output/test1.txt");
    //solarSystem.addFromFile("planetData.txt");

    //solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    //solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2.1 * M_PI, 0), 5.97219/ solarMass, "Earth");
    /*
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    solarSystem.addBody(vec3(-6.333487572394930E-02, -4.608453269808703E-01, -3.184761165078634E-02),
                        vec3(2.222816779156590E-02, -2.399853089908365E-03, -2.235205883702246E-03) * 365.25, 0.3302 / solarMass, "Mercury");

    vVerlet solver(dt);

    solarSystem.writeSystem();
    for (int i = 0; i < n; i++) {
        solarSystem.writePosVel();
        solver.step(solarSystem, 2);
        //solver.relativisticStep(solarSystem, 2);
    }
    */
    //solarSystem.recenter();
    //solarSystem.nullifyMomentum();

    //Euler solver(5.0 / n);

    // Mercury Orbit with Relativistic Gravity
    n = 100000000;
    totalTime = 100;
    dt = totalTime / n;
    SolarSystem mercurySystem("Output/mercury.txt");
    mercurySystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    mercurySystem.addBody(vec3(0.3075, 0, 0),
                          vec3(0, 12.44, 0), 0.3302 / solarMass, "Mercury");
    vVerlet solver(dt);
    mercurySystem.writeSystem();
    int i = 0;
    while (i < 250000) {
        mercurySystem.writePos();
        solver.relativisticStep(mercurySystem, 2);
        i++;
    }
    while (i < n - 250000) {
        solver.relativisticStep(mercurySystem, 2);
        i++;
    }
    while (i < n) {
        mercurySystem.writePos();
        solver.relativisticStep(mercurySystem, 2);
        i++;
    }

    return 0;
}
