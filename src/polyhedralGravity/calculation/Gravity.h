#pragma once

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityResult.h"
#include "polyhedralGravity/util/UtilityConstants.h"
#include "polyhedralGravity/util/UtilityContainer.h"
#include "spdlog/spdlog.h"

namespace polyhedralGravity {

    class Gravity {

        /**
         * The Polyhedron for which we evaluated the gravity at point P
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

        /**
         * Executes the whole calculation.
         * TODO: return values or via getResult()?
         */
        void calculate();

        /**
         * Calculates the G_ij vectors according to Tsoulis equation (18). These vectors are required to further
         * compute the plane unit and segment unit normals.
         *
         * The dimension of i will be equal to the number of faces, whereas the dimension j will be equal to 3 as the
         * given polyhedral's faces always consist of three segments/ nodes (triangles).
         * @return G vectors
         */
        std::vector<std::array<std::array<double, 3>, 3>> calculateGij();

        /**
         * Calculate the N_i vectors according to Tsoulis equation (19).
         *
         * The dimension of i will be equal to the number of faces.
         * @return plane unit normals
         */
        std::vector<std::array<double, 3>> calculatePlaneUnitNormals();

        /**
         * Calculates the segment unit normals according to Tsoulis equation (20).
         *
         * The dimension of i will be equal to the number of faces, whereas the dimension j mirrors the number of
         * segments forming one face. Since we always use triangles, j will be 3.
         * @return segment unit normals
         */
        std::vector<std::array<std::array<double, 3>, 3>> calculateSegmentUnitNormals();


    };

}