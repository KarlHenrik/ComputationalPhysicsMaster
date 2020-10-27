#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "solarSystem.h"
#include "celestialBody.h"
#include "vVerlet.h"
using namespace std;

int main(int argc, char const *argv[]) {
    int n;
    double totalTime, dt;
    double solarMass = 1988500;
    SolarSystem solarSystem("Output/testing.txt");
    vVerlet solver(dt);

    // ------------------- Some planets that will be reused ---------------------------------------
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ solarMass, "Earth");
    CelestialBody earthEll(vec3(1, 0, 0), vec3(0, 5, 0), 5.97219/ solarMass, "Earth");
    CelestialBody earthReal(vec3(-1.663457613546699E-01, 9.691203921746886E-01, -4.125583172010008E-05),
                            vec3(-1.723919408870981E-02, -2.981520896064708E-03, 4.254600200473125E-07) * 365.25, 5.97219/ solarMass, "Earth");
    CelestialBody sun(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun");
    CelestialBody jupyter(vec3(5.261470562232079E-01, -5.201022508399864E+00, 9.830503793253315E-03),
                          vec3(7.423674451944973E-03, 1.116865602956755E-03, -1.707572410053752E-04) * 365.25, 1898.13 / solarMass, "Jupyter");
    CelestialBody mercury(vec3(-6.333487572394930E-02, -4.608453269808703E-01, -3.184761165078634E-02),
                          vec3(2.222816779156590E-02, -2.399853089908365E-03, -2.235205883702246E-03) * 365.25, 0.3302 / solarMass, "Mercury");

    // ------------------- Cicrular Earth Sun with different exponents -----------------------------
    n = 100000;
    totalTime = 10;
    dt = totalTime / n;

    vector<double> exponents1{2, 2.1, 2.5, 3};
    for (double exponent : exponents1) {
        string ofile = "Output/Circular_exponent" + to_string(exponent) + ".txt";
        solarSystem.setofile(ofile);
        solarSystem.addBodyObj(earthCirc);
        solarSystem.addBodyObj(sun, true);
        solver.setdt(dt);
        solarSystem.writeSystem();
        for (int i = 0; i < n; i++) {
            solarSystem.writePosVel();
            solver.step(solarSystem, exponent);
        }
    }

    // ------------------- Elliptical Earth Sun with different exponents -----------------------------
    n = 100000;
    totalTime = 2;
    dt = totalTime / n;

    vector<double> exponents2{2, 2.1, 2.5, 3};
    for (double exponent : exponents2) {
        string ofile = "Output/Elliptical_exponent" + to_string(exponent) + ".txt";
        solarSystem.setofile(ofile);
        solarSystem.addBodyObj(earthEll);
        solarSystem.addBodyObj(sun, true);
        solver.setdt(dt);
        solarSystem.writeSystem();
        for (int i = 0; i < n; i++) {
            solarSystem.writePosVel();
            solver.step(solarSystem, exponent);
        }
    }


    // ------------------- Earth Escape Velocity -----------------------------
    n = 200000;
    totalTime = 100.0;
    dt = totalTime / n;

    vector<double> velocities{8, 8.8, 8.9, 9};

    for (double vel : velocities) {
        string ofile = "Output/Escape_velocity" + to_string(vel) + ".txt";
        solarSystem.setofile(ofile);
        solarSystem.addBody(vec3(1, 0, 0), vec3(vel, 0, 0), 5.97219 / solarMass, "Earth");
        solarSystem.addBodyObj(sun, true);
        solver.setdt(dt);
        solarSystem.writeSystem();
        for (int i = 0; i < n; i++) {
            solarSystem.writePos();
            solver.step(solarSystem);
        }
    }

    // ------------------- Only Earth -----------------------------
    n = 1000000;
    totalTime = 100;
    dt = totalTime / n;

    string ofile = "Output/onlyearth.txt";
    solarSystem.setofile(ofile);
    solarSystem.addBodyObj(earthReal);
    solarSystem.addBodyObj(sun, true);
    solver.setdt(dt);
    solarSystem.writeSystem();
    for (int i = 0; i < n; i++) {
        solarSystem.writePos();
        solver.step(solarSystem);

    }

    // ------------------- Earth and Jupyter -----------------------------
    n = 100000;
    totalTime = 10.0;
    dt = totalTime / n;

    vector<double> jmasses{1, 10, 1000};
    for (double jmass : jmasses) {
        string ofile = "Output/jupyterMass" + to_string(jmass) + ".txt";
        solarSystem.setofile(ofile);
        solarSystem.addBodyObj(earthReal);
        solarSystem.addBody(vec3(5.261470562232079E-01, -5.201022508399864E+00, 9.830503793253315E-03),
                            vec3(7.423674451944973E-03, 1.116865602956755E-03, -1.707572410053752E-04) * 365.25, 1898.13 * jmass / solarMass, "Jupyter");
        solarSystem.addBodyObj(sun, true);
        solver.setdt(dt);
        solarSystem.writeSystem();
        for (int i = 0; i < n; i++) {
            solarSystem.writePos();
            solver.step(solarSystem);
        }
    }

    // ------------------- Earth and Jupyter moving sun -----------------------------
    n = 100000;
    totalTime = 10;
    dt = totalTime / n;

    vector<double> jmasses2{1, 10, 1000};
    for (double jmass : jmasses2) {
        string ofile = "Output/jupyterMassMove" + to_string(jmass) + ".txt";
        solarSystem.setofile(ofile);
        solarSystem.addBodyObj(earthReal);
        solarSystem.addBody(vec3(5.261470562232079E-01, -5.201022508399864E+00, 9.830503793253315E-03),
                            vec3(7.423674451944973E-03, 1.116865602956755E-03, -1.707572410053752E-04) * 365.25, 1898.13 * jmass / solarMass, "Jupyter");
        solarSystem.addBodyObj(sun);
        solver.setdt(dt);
        solarSystem.writeSystem();
        for (int i = 0; i < n; i++) {
            solarSystem.writePos();
            solver.step(solarSystem);
        }
    }


    // ------------------ All planets long simulation --------------------
    n = 1000000;
    totalTime = 100;
    dt = totalTime / n;

    solarSystem.setofile("Output/allplong.txt");
    solarSystem.addFromFile("planetData.txt", 9);

    solver.setdt(dt);
    solarSystem.recenter();
    solarSystem.nullifyMomentum();
    solarSystem.writeSystem();
    for (int i = 0; i < n; i++) {
        solarSystem.writePos();
        solver.step(solarSystem);
    }


    // --------------Mercury Orbit with Relativistic Gravity----------------
    n = 7.5E7;
    totalTime = 100;
    dt = totalTime / n;

    solarSystem.setofile("Output/mercury.txt");
    solarSystem.addBodyObj(mercury);
    solarSystem.addBodyObj(sun, true);

    solver.setdt(dt);
    solarSystem.writeSystem();
    int i = 0;
    int inds = n / 400;
    while (i < inds) {
        solarSystem.writePos();
        solver.relativisticStep(solarSystem);
        i++;
    }
    while (i < n - inds) {
        solver.relativisticStep(solarSystem);
        i++;
    }
    while (i < n) {
        solarSystem.writePos();
        solver.relativisticStep(solarSystem);
        i++;
    }

    // --------------Mercury Orbit with Normal Gravity----------------
    n = 7.5E7;
    totalTime = 100;
    dt = totalTime / n;

    solarSystem.setofile("Output/mercuryNormal.txt");
    solarSystem.addBodyObj(mercury);
    solarSystem.addBodyObj(sun, true);

    solver.setdt(dt);
    solarSystem.writeSystem();
    int i2 = 0;
    int inds2 = n / 400;
    while (i2 < inds2) {
        solarSystem.writePos();
        solver.step(solarSystem);
        i2++;
    }
    while (i2 < n - inds2) {
        solver.step(solarSystem);
        i2++;
    }
    while (i2 < n) {
        solarSystem.writePos();
        solver.step(solarSystem);
        i2++;
    }

    return 0;
}
