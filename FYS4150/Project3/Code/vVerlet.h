#ifndef vVerlet_H
#define vVerlet_H

#include <iostream>
#include "vec3.h"
#include "solarSystem.h"

class vVerlet {
    public:
        double dt;
        vVerlet(double dt);
        void step(SolarSystem &ssys);
};

#endif // vVerlet_H
