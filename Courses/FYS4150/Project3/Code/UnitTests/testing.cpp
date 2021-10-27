#include "catch.hpp"

#include "../solarSystem.h"
#include "../vVerlet.h"
#include "../celestialBody.h"
#include "../euler.h"

TEST_CASE("Lager solarSystem objekt og legger til planeter", "[addBody, addBodyObj]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc, true);
    REQUIRE( solarSystem.bodies[1].vel[1] == 2 * M_PI );
    REQUIRE( solarSystem.bodies[1].pos[0] == 1 );
}

TEST_CASE("Oppdaterer posisjoner og hastigheter i solsystem med vVerlet", "[step]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc);

    int n = 1000; // antall tidssteg
    double totalTime = 1; // total tid
    double dt = totalTime / n; // størrelse på tidssteg
    vVerlet solver(dt); // vVerlet objekt med tidssteg dt
    for (int i = 0; i < n; i++) {
        solver.step(solarSystem); // finner akselerasjoner og oppdaterer posisjoner og hastigheter
    }
    REQUIRE( solarSystem.bodies[1].vel[1] != 2 * M_PI );
    REQUIRE( solarSystem.bodies[1].pos[0] != 1 );
}

TEST_CASE("Oppdaterer posisjoner og hastigheter i solsystem med euler", "[step]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc);

    int n = 1000; // antall tidssteg
    double totalTime = 1; // total tid
    double dt = totalTime / n; // størrelse på tidssteg
    Euler solver(dt); // vVerlet objekt med tidssteg dt
    for (int i = 0; i < n; i++) {
        solver.step(solarSystem); // finner akselerasjoner og oppdaterer posisjoner og hastigheter
    }
    REQUIRE( solarSystem.bodies[1].vel[1] != 2 * M_PI );
    REQUIRE( solarSystem.bodies[1].pos[0] != 1 );
}

TEST_CASE("tester effekt av nullifyMomentum med ingen parametre", "[nullifyMomentum]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc, true);

    vec3 vel0 = solarSystem.bodies[0].vel;
    solarSystem.nullifyMomentum();
    vec3 totalMomentum = vec3(0, 0, 0);
    for (CelestialBody &body : solarSystem.bodies) {
        totalMomentum += body.vel * body.mass;
    }
    REQUIRE( vel0[1] != solarSystem.bodies[0].vel[1] ); //y velocity of sun should change
    REQUIRE( totalMomentum[0] == Approx(0).epsilon(0.001) );
    REQUIRE( totalMomentum[1] == Approx(0).epsilon(0.001) );
    REQUIRE( totalMomentum[2] == Approx(0).epsilon(0.001) );
}

TEST_CASE("tester effekt av nullifyMomentum med parameter -1", "[nullifyMomentum]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc, true);
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen

    vec3 vel0 = solarSystem.bodies[3].vel;
    solarSystem.nullifyMomentum(-1);

    vec3 totalMomentum = vec3(0, 0, 0);
    for (CelestialBody &body : solarSystem.bodies) {
        totalMomentum += body.vel * body.mass;
    }
    REQUIRE( vel0[1] != solarSystem.bodies[3].vel[1] ); //y velocity of sun should change
}

TEST_CASE("tester effekt recenter", "[recenter]") {
    SolarSystem solarSystem("../Output/testing.txt"); // SolarSystem objekt med output fil
    solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
    CelestialBody earthCirc(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth");
    solarSystem.addBodyObj(earthCirc, true);
    solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen

    solarSystem.recenter();

    double totalMass = 0;
    for (CelestialBody &body : solarSystem.bodies) {
        totalMass += body.mass;
    }
    vec3 com = vec3(0, 0, 0); // center of mass
    for (CelestialBody &body : solarSystem.bodies) {
        com += body.pos * body.mass / totalMass;
    }
    REQUIRE( com[0] == Approx(0).epsilon(0.001) );
    REQUIRE( com[1] == Approx(0).epsilon(0.001) );
    REQUIRE( com[2] == Approx(0).epsilon(0.001) );
}
