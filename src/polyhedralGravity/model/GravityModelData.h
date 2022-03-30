#pragma once

#include <array>
#include <ostream>

namespace polyhedralGravity {

    /**
     * Contains the 3D distances l1_pq and l2_pq between P and the endpoints of segment pq and
     * the 1D distances s1_pq and s2_pq between P'' and the segment endpoints.
     * @note This struct is basically a named tuple
     */
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

        friend std::ostream &operator<<(std::ostream &os, const Distance &distance) {
            os << "l1: " << distance.l1 << " l2: " << distance.l2 << " s1: " << distance.s1 << " s2: " << distance.s2;
            return os;
        }
    };

    /**
     * Contains the Transcendental Expressions LN_pq and AN_pq for a given line segment pq of the polyhedron.
     * @note This struct is basically a named tuple
     */
    struct TranscendentalExpression {
        double ln;
        double an;

        bool operator==(const TranscendentalExpression &rhs) const {
            return ln == rhs.ln &&
                   an == rhs.an;
        }

        bool operator!=(const TranscendentalExpression &rhs) const {
            return !(rhs == *this);
        }

        friend std::ostream &operator<<(std::ostream &os, const TranscendentalExpression &expression) {
            os << "ln: " << expression.ln << " an: " << expression.an;
            return os;
        }
    };


    /**
     * A struct describing a plane in Hessian Normal Form:
     * ax + by + cz + d = 0
     * where a,b,c are the plane's normal
     * and d as the signed distance to the plane from the origin along the normal.
     */
    struct HessianPlane {
        double a;
        double b;
        double c;
        double d;

        bool operator==(const HessianPlane &rhs) const {
            return a == rhs.a &&
                   b == rhs.b &&
                   c == rhs.c &&
                   d == rhs.d;
        }

        bool operator!=(const HessianPlane &rhs) const {
            return !(rhs == *this);
        }
    };

    /**
     * A data structure containing the result of the polyhedral gravity model's evaluation.
    */
    class GravityModelResult {

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
        GravityModelResult()
                : p{0, 0, 0} {}

        /**
         * Creates new empty Result for a specific point P
         * @param p1 - point p
         */
        explicit GravityModelResult(const std::array<double, 3> &p1)
                : p{p1} {}

    };

}
