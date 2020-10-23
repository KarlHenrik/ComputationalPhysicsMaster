#ifndef vVerlet_H
#define vVerlet_H

#include <iostream>
#include "vec3.h"
#include "solarSystem.h"

class vVerlet {
    public:
        double dt;
        double dt2;
        vVerlet(double dt);
        void step(SolarSystem &ssys, double exponent = 2);
        void relativisticStep(SolarSystem &ssys, double exponent = 2);
};

#endif // vVerlet_H
