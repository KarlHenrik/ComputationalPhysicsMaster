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
}

TEST_CASE("Oppdaterer posisjoner i solsystem med vVerlet", "[step]") {
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
}

TEST_CASE("test case description", "[functionName]") {
    REQUIRE( 1 == 1 );
}
