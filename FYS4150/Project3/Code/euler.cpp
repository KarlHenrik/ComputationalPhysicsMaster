#include "euler.h"

Euler::Euler(double dt_) {
    dt = dt_;
}

void Euler::step(SolarSystem &ssys, double exponent) {
    ssys.calcAcc(exponent);

    for(CelestialBody &body : ssys.bodies) {
        body.pos += body.vel * dt;
        body.vel += body.acc * dt;
    }
}
