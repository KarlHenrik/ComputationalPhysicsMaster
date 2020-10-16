#include "solarSystem.h"

SolarSystem::SolarSystem(string ofilename) {
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
}

void SolarSystem::setofile(string ofilename) {
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
}

void SolarSystem::addBody(vec3 pos, vec3 vel, double mass, string name) {
    double solarMass = 1988500;
    bodies.push_back( CelestialBody(pos, vel, mass/solarMass, name) );
    nbodies ++;
}

void SolarSystem::addFromFile(string infilename, int nread) {
    ifstream infile;
    infile.open(infilename);
    double mass, x, y, z, vx, vy, vz;
    string name, line;
    vector<string> splitLine;
    string temp1, temp2, temp3;
    int nadded = 0;
    while (getline(infile, line))
    {
        splitLine = split(line);
        name = splitLine[0];
        mass = sci2dbl(splitLine[1]);

        getline(infile, line);
        splitLine = split(line);
        x = sci2dbl(splitLine[0]);
        y = sci2dbl(splitLine[1]);
        z = sci2dbl(splitLine[2]);

        getline(infile, line);
        splitLine = split(line);
        vx = sci2dbl(splitLine[0]);
        vy = sci2dbl(splitLine[1]);
        vz = sci2dbl(splitLine[2]);

        addBody(vec3(x,y,z), vec3(vx,vy,vz) * 365.25, mass, name);
        nadded++;
        if(nadded == nread) {
            return;
        }
    }
}

void SolarSystem::calcAcc(double exponent) {
    // Set all accelerations to 0
    for (int i = 0; i < nbodies; i++) {
        bodies[i].acc = vec3(0, 0, 0);
    }
    double g = 4 * M_PI * M_PI; // gravitational constant for units AU, year and solar mass
    vec3 pos1to2;
    double r2, forceFactor;
    for (int i = 0; i < nbodies; i++) {
        CelestialBody &body1 = bodies[i];
        for (int j = i + 1; j < nbodies; j++) { // we only look at each pair of planets once
            CelestialBody &body2 = bodies[j];
            pos1to2 = body2.pos - body1.pos;
            r2 = pos1to2.lengthSquared();
            forceFactor = g / (r2 * sqrt(r2));
            body1.acc += forceFactor * body2.mass * pos1to2;
            body2.acc += -forceFactor * body1.mass * pos1to2;
        }
    }
}

void SolarSystem::writeSystem() {
    return;
}
void SolarSystem::writePosVel() {
    for (CelestialBody &body : bodies) {
        ofile << setw(15) << setprecision(8) << body.pos[0];
        ofile << setw(15) << setprecision(8) << body.pos[1];
        ofile << setw(15) << setprecision(8) << body.pos[2];
    }
    ofile << endl;
    return;
}

void SolarSystem::recenter() {
    return;
}
void SolarSystem::nullifyMomentum(int bodyIndex) {
    return;
}

double SolarSystem::sci2dbl(string s) {
    istringstream iss(s);
    double d;
    iss >> d;
    return d;
}

vector<string> SolarSystem::split(string line) {
    vector<string> splitLine;
    istringstream iss(line);
    for(string line; iss >> line; )
        splitLine.push_back(line);
    return splitLine;
}
