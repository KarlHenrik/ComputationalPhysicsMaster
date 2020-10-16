#include "celestialBody.h"

CelestialBody::CelestialBody(vec3 pos_, vec3 vel_, double mass_, string name_) {
    pos = pos_;
    vel = vel_;
    mass = mass_;
    name = name_;
}
