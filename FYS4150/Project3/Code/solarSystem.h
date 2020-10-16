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

using namespace std;

class SolarSystem {
    public:
        int nbodies = 0;
        vector<CelestialBody> bodies;
        ofstream ofile;

        SolarSystem(string ofilename);
        void setofile(string ofilename);
        void addBody(vec3 pos, vec3 vel, double mass, string name);
        void addFromFile(string infilename, int nread = 10);

        void calcAcc(double exponent = 2);

        void writeSystem();
        void writePosVel();

        void recenter();
        void nullifyMomentum(int bodyIndex = 0);

        double sci2dbl(string s);
        vector<string> split(string line);
};

#endif // SolarSystem_H
