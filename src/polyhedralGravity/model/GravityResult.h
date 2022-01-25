#pragma once

#include <array>

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

    GravityResult()
            : p{0, 0, 0} {}

    explicit GravityResult(const std::array<double, 3> &p1)
            : p{p1} {}

};