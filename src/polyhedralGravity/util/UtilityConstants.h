#pragma once

namespace polyhedralGravity::util {

    /**
     * The epsilon used in the polyhedral gravity model.
     * @related Used for the sgn() function to determine the sign of a double value. Different compilers
     * produce different results if no epsilon is applied for the comparison!
     */
    constexpr double epsilon = 1e-15;

    /**
     * The gravitational constant G in [m^3/(kg*s^2)].
     * @related in his paper above Equation (4)
     */
    constexpr double gravitationalConstant = 6.67430e-11;

    /**
     * The assumed constant density rho for a polyhedron after Tsoulis paper in [kg/m^3].
     * @related in his paper above Equation (4)
     */
    constexpr double defaultConstantDensity = 2670.0;

}
