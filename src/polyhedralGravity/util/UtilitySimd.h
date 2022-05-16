#pragma once

#include <array>
#include <utility>

#include "xsimd/xsimd.hpp"

namespace polyhedralGravity::util {

    /**
     * Applies the euclidean norm/ L2-norm to a std::array of size 3 in a efficient manner.
     * @param array - array of size 3
     * @return an double containing the L2 norm
     */
    inline double euclideanNorm(const std::array<double, 3> array) {
        xsimd::batch<double> reg1({array[0], array[1]});
        reg1 = xsimd::mul(reg1, reg1);
        double reg3 = array[2] * array[2];
        return std::sqrt(reg1.get(0) + reg1.get(1) + reg3);
    }

    inline std::pair<double, double> euclideanNorm(
            const std::array<double, 3> array1, const std::array<double, 3> array2) {
        xsimd::batch<double> reg1({array1[0], array1[1]});
        xsimd::batch<double> reg2({array1[2], array2[0]});
        xsimd::batch<double> reg3({array2[1], array2[2]});
        reg1 = xsimd::mul(reg1, reg1);
        reg2 = xsimd::mul(reg2, reg2);
        reg2 = xsimd::mul(reg3, reg3);

        reg1 = xsimd::batch<double>({reg1.get(0) + reg1.get(1) + reg2.get(0), reg2.get(0) + reg3.get(0) + reg3.get(1)});
        reg1 = xsimd::sqrt(reg1);
        return std::make_pair(reg1.get(0), reg1.get(1));
    }

    inline std::array<double, 3> euclideanNorm(const std::array<double, 3> array1,
                                               const std::array<double, 3> array2, const std::array<double, 3> array3) {
        auto pair = euclideanNorm(array2, array3);
        return {euclideanNorm(array1), pair.first, pair.second};
    }

    inline std::array<double, 3> operator-(const std::array<double, 3> array1, const std::array<double, 3> array2) {
        xsimd::batch<double> reg1({array1[0], array1[1]});
        xsimd::batch<double> reg2({array2[0], array2[1]});
        reg1 = xsimd::sub(reg1, reg2);
        double reg3 = array1[2] - array2[2];
        return {reg1.get(0), reg1.get(1), reg3};
    }

    inline std::array<double, 3> operator+(const std::array<double, 3> array1, const std::array<double, 3> array2) {
        xsimd::batch<double> reg1({array1[0], array1[1]});
        xsimd::batch<double> reg2({array2[0], array2[1]});
        reg1 = xsimd::add(reg1, reg2);
        double reg3 = array1[2] + array2[2];
        return {reg1.get(0), reg1.get(1), reg3};
    }

}