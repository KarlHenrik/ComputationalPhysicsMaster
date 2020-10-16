#ifndef CelestialBody_H
#define CelestialBody_H

#include "vec3.h"

using namespace std;

class CelestialBody {
    public:
        vec3 pos;
        vec3 vel;
        vec3 acc;
        double mass;
        string name;

        CelestialBody(vec3 pos, vec3 vel, double mass, string name);
};

#endif // CelestialBody_H
