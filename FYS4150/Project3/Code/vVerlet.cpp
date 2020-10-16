#include "vVerlet.h"

vVerlet::vVerlet(double dt_) {
    dt = dt_;
}

void vVerlet::step(SolarSystem &ssys) {
    ssys.calcAcc();

    for(CelestialBody &body : ssys.bodies) {
        body.vel += body.acc * dt;
        body.pos += body.vel * dt;
    }
}
