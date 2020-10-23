#include "vVerlet.h"

vVerlet::vVerlet(double dt_) {
    dt = dt_;
    dt2 = dt * dt;
}

void vVerlet::step(SolarSystem &ssys, double exponent) {
    // update using old acceleration
    for(CelestialBody &body : ssys.bodies) {
        body.pos += dt * body.vel + dt2 / 2 * body.acc;
        body.vel += dt / 2 * body.acc;
    }
    ssys.calcAcc(exponent);
    // update velocity further using new accelatation
    for(CelestialBody &body : ssys.bodies) {
        body.vel += dt / 2 * body.acc;
    }
}

void vVerlet::relativisticStep(SolarSystem &ssys, double exponent) {
    CelestialBody &body = ssys.bodies[1];
    vec3 oldAcc = body.acc;
    body.pos += dt * body.vel + dt2 / 2 * body.acc;
    ssys.relativisticAcc(exponent);
    body.vel += dt / 2 * oldAcc;
    body.vel += dt / 2 * body.acc;
}
