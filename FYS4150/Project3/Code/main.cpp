#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarSystem.h"
#include "euler.h"
#include "vVerlet.h"
using namespace std;

int main(int argc, char const *argv[]) {
    /*
    // ------------------ All planets simulation --------------------
    int n = 100000;
    double totalTime = 10;
    double dt = totalTime / n;
    double solarMass = 1988500;

    SolarSystem solarSystem("Output/allp.txt");
    solarSystem.addFromFile("planetData.txt");

    vVerlet solver(dt);
    solarSystem.recenter();
    solarSystem.nullifyMomentum();
    solarSystem.writeSystem();
    for (int i = 0; i < n; i++) {
        solarSystem.writePosVel();
        solver.step(solarSystem, 2);
    }
    */

    /*
    // ------------------- Earth Sun -----------------------------
    int n = 100000;
    double totalTime = 10;
    double dt = totalTime / n;
    double solarMass = 1988500;

    SolarSystem solarSystem("Output/earthSun.txt");
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ solarMass, "Earth");
    vVerlet solver(dt);
    solarSystem.writeSystem();
    for (int i = 0; i < n; i++) {
        solarSystem.writePosVel();
        solver.step(solarSystem, 2);
    }
    */



    // --------------Mercury Orbit with Relativistic Gravity----------------
    int n = 7.5E7;
    double totalTime = 100;
    double dt = totalTime / n;
    double solarMass = 1988500;

    SolarSystem mercurySystem("Output/mercury.txt");
    mercurySystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    //mercurySystem.addBody(vec3(0.3075, 0, 0),
    //                      vec3(0, 12.44, 0), 0.3302 / solarMass, "Mercury");
    mercurySystem.addBody(vec3(-6.333487572394930E-02, -4.608453269808703E-01, -3.184761165078634E-02),
                        vec3(2.222816779156590E-02, -2.399853089908365E-03, -2.235205883702246E-03) * 365.25, 0.3302 / solarMass, "Mercury");

    vVerlet solver(dt);
    mercurySystem.writeSystem();
    int i = 0;
    int inds = n / 400;
    while (i < inds) {
        mercurySystem.writePos();
        solver.relativisticStep(mercurySystem, 2);
        i++;
    }
    while (i < n - inds) {
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
