#include "euler.h"

Euler::Euler(double dt_) {
    dt = dt_;
}

void Euler::setdt(double dt_) {
    dt = dt_;
}

void Euler::step(SolarSystem &ssys, double exponent) {
    // update using old acceleration
    ssys.calcAcc(exponent);
    for (int i = 0; i < ssys.nbodiesMoveable; i++) {
        CelestialBody &body = ssys.bodies[i];
        body.pos += dt * body.vel;
        body.vel += dt * body.acc;
    }
}
