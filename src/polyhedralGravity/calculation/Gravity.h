#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityResult.h"
#include "spdlog/spdlog.h"

class Gravity {

    /**
     * The Polyhedron for which to evaluated the gravity at point P
     */
    const Polyhedron _polyhedron;

    /**
     * The constant density of the polyhedron in [kg/m^3].
     */
    const double _density;

    /**
     * The result of the evaluation of the gravity model
     */
    GravityResult _gravityResult{};

public:

    /**
     * Construct as new Gravity Calculation with an polyhedron as input
     * @param polyhedron - Poylhedron
     * The density is initialized with the default density for asteroids (2 g/cm^3 = 2000 kg/m^3).
     */
    explicit Gravity(const Polyhedron &polyhedron)
            : _polyhedron{polyhedron},
              _density{2000.0} {}

    /**
     * Construct as new Gravity Calculation with an polyhedron as input
     * @param polyhedron - Poylhedron
     * @param density - the constant density of the Polyhedron
     */
    Gravity(const Polyhedron &polyhedron, double density)
            : _polyhedron{polyhedron},
              _density{density} {}


    void calculate();


};
