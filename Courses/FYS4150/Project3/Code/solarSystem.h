#ifndef SolarSystem_H
#define SolarSystem_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "vec3.h"
#include "celestialBody.h"
#include <iostream>

using namespace std;

class SolarSystem {
    public:
        double g;
        int nbodies = 0;
        int nbodiesMoveable = 0;
        vector<CelestialBody> bodies;
        ofstream ofile;

        SolarSystem(string ofilename);
        void setofile(string ofilename);
        void addBody(vec3 pos, vec3 vel, double mass, string name, bool ignore = false);
        void addBodyObj(CelestialBody body, bool ignore = false);
        void addFromFile(string infilename, int nread = 10);

        void calcAcc(double exponent = 2);
        void relativisticAcc(double exponent = 2);

        void writeSystem();
        void writePos();
        void writePosVel();

        void recenter();
        void nullifyMomentum(int bodyIndex = 0);

        double sci2dbl(string s);
        vector<string> split(string line);
};

#endif // SolarSystem_H
