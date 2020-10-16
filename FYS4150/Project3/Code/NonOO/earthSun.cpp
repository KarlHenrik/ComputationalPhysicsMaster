#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "time.h"

#include "../vec3.h"

using namespace std;

double h;

void solveTime(int);
void solveWrite(string, int);
void accel(vec3&, vec3);
void euler(vec3&, vec3&, vec3&);
void vverlet(vec3&, vec3&, vec3&);

int main(int argc, char const *argv[]) {
    for (int i = 2; i < 6; i++) {
        solveWrite("euler", i);
        solveWrite("vverlet", i);
    }
    solveTime(8);
    return 0;
}

void solveTime(int imax) {
    // Setting up output file
    string ofilename = "Output/earthSun/earthSunTiming.txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << "n";
    ofile << setw(15) << "Velocity Verlet";
    ofile << setw(15) << "Euler" << endl;

    clock_t start1, finish1, start2, finish2;
    // Solving equations
    for (int i = 5; i <= imax; i++) {
        // Parameters for solver
        h = pow(10, -i);
        int n = 5 / h;
        // Timing Velocity Verlet
        vec3 earthPos(1, 0, 0);
        vec3 earthVel(0, 2 * M_PI, 0);
        vec3 earthAcc(0, 0, 0);
        start1 = clock();
        for (int i = 0; i < n; i++) {
            vverlet(earthAcc, earthVel, earthPos);
        }
        finish1 = clock();

        // Timing Euler
        vec3 earthPos2(1, 0, 0);
        vec3 earthVel2(0, 2 * M_PI, 0);
        vec3 earthAcc2(0, 0, 0);
        start2 = clock();
        for (int i = 0; i < n; i++) {
            euler(earthAcc2, earthVel2, earthPos2);
        }
        finish2 = clock();

        ofile << setw(15) << setprecision(8) << n;
        ofile << setw(15) << setprecision(8) << (finish1 - start1) / ((double) CLOCKS_PER_SEC);
        ofile << setw(15) << setprecision(8) << (finish2 - start2) / ((double) CLOCKS_PER_SEC) << endl;
    }
}


void solveWrite(string solver, int i) {
    // Setting up output file
    string ofilename = "Output/earthSun/earthSun_" + solver + to_string(i) + ".txt";
    ofstream ofile;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    // i     x       y      <
    ofile << setw(15) << "x";
    ofile << setw(15) << "y";
    ofile << setw(15) << "z";
    ofile << setw(15) << "vx";
    ofile << setw(15) << "vy";
    ofile << setw(15) << "vz" << endl;
    // Initial conditions
    vec3 earthPos(1, 0, 0);
    vec3 earthVel(0, 2 * M_PI, 0);
    vec3 earthAcc(0, 0, 0);

    // Parameters for solver
    h = pow(10, -i);
    int n = 5 / h;

    // Solving equations
    for (int i = 0; i < n; i++) {
        // Calculation next iteration
        if (solver == "euler") {
            euler(earthAcc, earthVel, earthPos);
        } else {
            vverlet(earthAcc, earthVel, earthPos);
        }
        // Writing position and velocity to file
        ofile << setw(15) << setprecision(8) << earthPos[0];
        ofile << setw(15) << setprecision(8) << earthPos[1];
        ofile << setw(15) << setprecision(8) << earthPos[2];
        ofile << setw(15) << setprecision(8) << earthVel[0];
        ofile << setw(15) << setprecision(8) << earthVel[1];
        ofile << setw(15) << setprecision(8) << earthVel[2] << endl;
    }
}

void accel(vec3 &earthAcc, vec3 earthPos) {
    double g = 4 * M_PI * M_PI;
    double r2 = earthPos.lengthSquared();
    earthAcc = -g / (r2 * sqrt(r2)) * earthPos;
}

void euler(vec3 &earthAcc, vec3 &earthVel, vec3 &earthPos) {
    accel(earthAcc, earthPos);
    earthPos += h * earthVel;
    earthVel += h * earthAcc;
}

void vverlet(vec3 &earthAcc, vec3 &earthVel, vec3 &earthPos) {
    earthPos += h * earthVel + h * h / 2 * earthAcc;
    earthVel += h / 2 * earthAcc;
    accel(earthAcc, earthPos);
    earthVel += h / 2 * earthAcc;
}
