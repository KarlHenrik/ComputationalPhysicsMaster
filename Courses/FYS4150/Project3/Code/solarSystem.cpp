#include "solarSystem.h"

SolarSystem::SolarSystem(string ofilename) {
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
}

void SolarSystem::setofile(string ofilename) {
    ofile.close();
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    nbodies = 0;
    nbodiesMoveable = 0;
    bodies.clear();
}

void SolarSystem::addBody(vec3 pos, vec3 vel, double mass, string name, bool ignore) {
    bodies.push_back( CelestialBody(pos, vel, mass, name) );
    nbodies++;
    if (ignore == false) { // the ignored bodies must be added last, and will not be moved or have their positions printed
        nbodiesMoveable++;
    }
}

void SolarSystem::addBodyObj(CelestialBody body, bool ignore) {
    bodies.push_back( body );
    nbodies++;
    if (ignore == false) { // the ignored bodies must be added last, and will not be moved or have their positions printed
        nbodiesMoveable++;
    }
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

        addBody(vec3(x,y,z), vec3(vx,vy,vz) * 365.25, mass / 1988500.0, name);
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
    double r, forceFactor, rExp;
    for (int i = 0; i < nbodies; i++) {
        CelestialBody &body1 = bodies[i];
        for (int j = i + 1; j < nbodies; j++) { // we only look at each pair of planets once
            CelestialBody &body2 = bodies[j];
            pos1to2 = body2.pos - body1.pos;
            r = pos1to2.length();
            rExp = pow(r, exponent);
            forceFactor = g / (rExp * r);
            body1.acc += forceFactor * body2.mass * pos1to2;
            body2.acc += -forceFactor * body1.mass * pos1to2;
        }
    }
}

void SolarSystem::writeSystem() {
    for (CelestialBody &body : bodies) {
        ofile << setw(15) << setprecision(8) << body.name;
        ofile << setw(15) << setprecision(8) << body.mass;
    }
    ofile << endl;
    return;
}

void SolarSystem::writePos() {
    for (int i = 0; i < nbodiesMoveable; i++) {
        CelestialBody &body = bodies[i];
        ofile << setw(15) << setprecision(8) << body.pos[0];
        ofile << setw(15) << setprecision(8) << body.pos[1];
        ofile << setw(15) << setprecision(8) << body.pos[2];
    }
    ofile << endl;
    return;
}

void SolarSystem::writePosVel() {
    for (int i = 0; i < nbodiesMoveable; i++) {
        CelestialBody &body = bodies[i];
        ofile << setw(15) << setprecision(8) << body.pos[0];
        ofile << setw(15) << setprecision(8) << body.pos[1];
        ofile << setw(15) << setprecision(8) << body.pos[2];

        ofile << setw(15) << setprecision(8) << body.vel[0];
        ofile << setw(15) << setprecision(8) << body.vel[1];
        ofile << setw(15) << setprecision(8) << body.vel[2];
    }
    ofile << endl;
    return;
}

void SolarSystem::recenter() { //making the center of mass (0, 0, 0)
    double totalMass = 0;
    for (CelestialBody &body : bodies) {
        totalMass += body.mass;
    }
    vec3 com = vec3(0, 0, 0); // center of mass
    for (CelestialBody &body : bodies) {
        com += body.pos * body.mass / totalMass;
    }
    for (CelestialBody &body : bodies) {
        body.pos -= com;
    }
    return;
}

void SolarSystem::nullifyMomentum(int bodyIndex) { //making the total momentum 0
    vec3 totalMomentum = vec3(0, 0, 0);
    for (CelestialBody &body : bodies) {
        totalMomentum += body.vel * body.mass;
    }
    if (bodyIndex == -1) {
        bodies[nbodies].vel -= totalMomentum / bodies[nbodies].mass;
    }
    bodies[bodyIndex].vel -= totalMomentum / bodies[bodyIndex].mass; // changing the velicity of bodies[bodyIndex] so that the total momentum of the solar system is 0
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

void SolarSystem::relativisticAcc(double exponent) {
    // Set all accelerations to 0
    for (int i = 0; i < nbodies; i++) {
        bodies[i].acc = vec3(0, 0, 0);
    }
    double g = 4 * M_PI * M_PI; // gravitational constant for units AU, year and solar mass
    double c2 = 3999262982.49891169; //c^2 for c with units AU/year
    double r, forceFactor, rExp;

    CelestialBody &body1 = bodies[0];
    CelestialBody &body2 = bodies[1];
    r = body1.pos.length();
    rExp = pow(r, exponent);
    forceFactor = g / (rExp * r);
    body1.acc += -forceFactor * body2.mass * body1.pos * (1.0 + (3.0 * body1.pos.cross(body1.vel).lengthSquared() ) / (r * r * c2) );
}
