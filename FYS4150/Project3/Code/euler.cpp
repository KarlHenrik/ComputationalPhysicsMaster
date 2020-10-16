#include "euler.h"

Euler::Euler(double dt_) {
    dt = dt_;
}

void Euler::step(SolarSystem &ssys) {
    ssys.calcAcc();

    for(CelestialBody &body : ssys.bodies) {
        body.pos += body.vel * dt;
        body.vel += body.acc * dt;
    }
}
