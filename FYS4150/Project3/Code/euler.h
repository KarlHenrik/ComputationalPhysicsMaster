#ifndef Euler_H
#define Euler_H

#include <iostream>
#include "vec3.h"
#include "solarSystem.h"

class Euler {
    public:
        double dt;
        Euler(double dt);
        void step(SolarSystem &ssys, double exponent = 2);
};

#endif // Euler_H
