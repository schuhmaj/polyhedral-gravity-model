#pragma once

namespace polyhedralGravity::util {

    /**
     * The gravitational constant G in [m^3/(kg*s^2)].
     * @related in his paper above Equation (4)
     */
    constexpr double gravitationalConstant = 6.67259e-11;

    /**
     * The assumed constant density /rho for a polyhedron after Tsoulis paper in [kg/m^3].
     * @related in his paper above Equation (4)
     */
    constexpr double defaultConstantDensity = 2670.0;

}
