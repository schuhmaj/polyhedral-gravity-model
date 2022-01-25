#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityResult.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "spdlog/spdlog.h"

class Gravity {

    /**
     * The Polyhedron for which to evaluated the gravity at point P
     */
    const Polyhedron _polyhedron;

    /**
     * The constant density of the polyhedron in [kg/m^3].
     * The density is initialized with the default constant density 2670.0 from Tsoulis Paper (above (4)).
     */
    const double _density{util::defaultConstantDensity};

    /**
     * The result of the evaluation of the gravity model
     */
    GravityResult _gravityResult{};

public:

    /**
     * Construct as new Gravity Calculation with an polyhedron as input
     * @param polyhedron - Poylhedron
     */
    explicit Gravity(const Polyhedron &polyhedron)
            : _polyhedron{polyhedron} {}

    /**
     * Construct as new Gravity Calculation with an polyhedron as input
     * @param polyhedron - Poylhedron
     * @param density - the constant density of the Polyhedron
     */
    Gravity(const Polyhedron &polyhedron, double density)
            : _polyhedron{polyhedron},
              _density{density} {}


    void calculate();

    /**
     * Calculates the LN_pq values based on Equation (14).
     * The subscript q is the polyhedral segment of one face p of the complete polyhedron.
     * @param p - the index of the polyhedral face
     * @param q - the polyhedral segment
     */
    double calculateLNpq(size_t p , size_t q);

    std::vector<std::array<std::array<double, 3>, 3>> calculateGij();


};
