#ifndef euler_H
#define euler_H

#include <iostream>
#include "vec3.h"
#include "solarSystem.h"

class Euler {
    public:
        double dt;
        Euler(double dt);
        void setdt(double dt);
        void step(SolarSystem &ssys, double exponent = 2);
};

#endif // euler_H
