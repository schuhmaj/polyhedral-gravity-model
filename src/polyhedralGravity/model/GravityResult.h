#pragma once

#include <array>

namespace polyhedralGravity {

    struct Distance {
        double l1;
        double l2;
        double s1;
        double s2;

        bool operator==(const Distance &rhs) const {
            return l1 == rhs.l1 &&
                   l2 == rhs.l2 &&
                   s1 == rhs.s1 &&
                   s2 == rhs.s2;
        }

        bool operator!=(const Distance &rhs) const {
            return !(rhs == *this);
        }
    };

    struct TranscendentalExpression {
        double ln;
        double an;
    };

/**
 * A data structure containing the result of the polyhedral gravity model's evaluation.
 */
    class GravityResult {

    public:

        /**
         * The point P at which the gravity model was evaluated.
         */
        const std::array<double, 3> p;

        /**
         * The gravitational potential in [m^2/s^2] <--> [J/kg] at point P.
         * @related Equation (1) and (11) of Tsoulis Paper, here referred as V
         */
        double gravitationalPotential{};

        /**
         * The first order derivatives of the gravitational potential in [m/s^2].
         * The array contains the derivatives depending on the coordinates x-y-z in this order.
         * @related Equation (2) and (12) of Tsoulis Paper, here referred as Vx, Vy, Vz
         */
        std::array<double, 3> gravitationalPotentialDerivative{};

        /**
         * The second order derivatives or also called gradiometric Tensor in [1/s^2].
         * The array contains the second order derivatives in the following order xx, yy, zz, xy, xz, yz.
         * @related Equation (3) and (13) of Tsoulis Paper, here referred as Vxx, Vyy, Vzz, Vxy, Vxz, Vyz
         */
        std::array<double, 6> gradiometricTensor{};

        /**
         * Creates a new empty Result for the origin in P(0, 0, 0)
         */
        GravityResult()
                : p{0, 0, 0} {}

        /**
         * Creates new empty Result for a specific point P
         * @param p1 - point p
         */
        explicit GravityResult(const std::array<double, 3> &p1)
                : p{p1} {}

    };

}
