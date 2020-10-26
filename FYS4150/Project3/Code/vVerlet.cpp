#include "vVerlet.h"

vVerlet::vVerlet(double dt_) {
    dt = dt_;
    dto2 = dt_/2;
    dt2o2 = dt * dt / 2;
}

void vVerlet::setdt(double dt_) {
    dt = dt_;
    dto2 = dt_/2;
    dt2o2 = dt_ * dt_ / 2;
}

void vVerlet::step(SolarSystem &ssys, double exponent) {
    // update using old acceleration
    for (int i = 0; i < ssys.nbodiesMoveable; i++) {
        CelestialBody &body = ssys.bodies[i];
        body.pos += dt * body.vel + dt2o2 * body.acc;
        body.vel += dto2 * body.acc;
    }
    ssys.calcAcc(exponent);
    // update velocity further using new accelatation
    for (int i = 0; i < ssys.nbodiesMoveable; i++) {
        CelestialBody &body = ssys.bodies[i];
        body.vel += dto2 * body.acc;
    }
}

void vVerlet::relativisticStep(SolarSystem &ssys, double exponent) {
    CelestialBody &body = ssys.bodies[0];
    vec3 oldAcc = body.acc;
    body.pos += dt * body.vel + dt2o2 * body.acc;
    ssys.relativisticAcc(exponent);
    body.vel += dto2 * oldAcc;
    body.vel += dto2 * body.acc;
}
