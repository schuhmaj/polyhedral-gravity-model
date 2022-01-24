#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "spdlog/spdlog.h"

class Gravity {

    /**
     * Point P were the gravity should be evaluated.
     * Currently does only support coordinates (0, 0, 0)
     */
    const std::array<double, 3> _p{0, 0, 0};

    /**
     * The Polyhedron for which to evaluated the gravity at point P
     */
    const Polyhedron _polyhedron;

public:

    /**
     * Construct as new Gravity Calculation with an polyhedron as input
     * @param polyhedron - Poylhedron
     */
    explicit Gravity(const Polyhedron &polyhedron)
            : _polyhedron{polyhedron} {}


    void calculate();


};
