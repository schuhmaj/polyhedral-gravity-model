#pragma once

#include <array>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>
#include <string>
#include <iostream>

#include "xsimd/xsimd.hpp"

namespace polyhedralGravity::util {

    /**
     * Alias for two-dimensional array with size M and N.
     * M is the major size.
     */
    template<typename T, size_t M, size_t N>
    using Matrix = std::array<std::array<T, N>, M>;

    /**
     * Applies the Operation Minus to two Containers piece by piece.
     * @example {1, 2, 3} - {1, 1, 1} = {0, 1, 2}
     * @tparam Container
     * @param lhs - minuend
     * @param rhs - subtrahend
     * @return the difference
     */
    template<typename Container>
    inline Container operator-(const Container &lhs, const Container &rhs) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::batch<T> rBatch = xsimd::load_unaligned(&rhs[i]);
            xsimd::store_aligned(&result[i], lBatch - rBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] - rhs[i];
        }
        return result;
    }

    /**
    * Applies the Operation Plus to two Containers piece by piece.
    * @example {1, 2, 3} + {1, 1, 1} = {2, 3, 4}
    * @tparam Container
    * @param lhs - addend
    * @param rhs - addend
    * @return the sum
    */
    template<typename Container>
    inline Container operator+(const Container &lhs, const Container &rhs) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::batch<T> rBatch = xsimd::load_unaligned(&rhs[i]);
            xsimd::store_aligned(&result[i], lBatch + rBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] + rhs[i];
        }
        return result;
    }

    /**
    * Applies the Operation * to two Containers piece by piece.
    * @example {1, 2, 3} * {2, 2, 2} = {2, 4, 6}
    * @tparam Container
    * @param lhs - multiplicand
    * @param rhs - multiplicand
    * @return the product
    */
    template<typename Container>
    inline Container operator*(const Container &lhs, const Container &rhs) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::batch<T> rBatch = xsimd::load_unaligned(&rhs[i]);
            xsimd::store_aligned(&result[i], lBatch * rBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] * rhs[i];
        }
        return result;
    }

    /**
    * Applies the Operation / to two Containers piece by piece.
    * @example {1, 2, 3} * {1, 2, 3} = {1, 1, 1}
    * @tparam Container
    * @param lhs - multiplicand
    * @param rhs - multiplicand
    * @return the product
    */
    template<typename Container>
    inline Container operator/(const Container &lhs, const Container &rhs) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::batch<T> rBatch = xsimd::load_unaligned(&rhs[i]);
            xsimd::store_aligned(&result[i], lBatch / rBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] / rhs[i];
        }
        return result;
    }

    /**
    * Applies the Operation + to a Container and a Scalar.
    * @example {1, 2, 3} + 2 = {3, 4, 5}
    * @tparam Container
    * @tparam Scalar
    * @param lhs - addend
    * @param scalar - addend
    * @return a Container
    */
    template<typename Container, typename Scalar>
    inline Container operator+(const Container &lhs, const Scalar &scalar) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        xsimd::batch<T> scalarBatch = xsimd::broadcast(scalar);
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::store_aligned(&result[i], lBatch + scalarBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] + scalar;
        }
        return result;
    }

    /**
    * Applies the Operation - to a Container and a Scalar.
    * @example {1, 2, 3} - 2 = {-1, 0, 1}
    * @tparam Container
    * @tparam Scalar
    * @param lhs - minuend
    * @param scalar - subtrahend
    * @return a Container
    */
    template<typename Container, typename Scalar>
    inline Container operator-(const Container &lhs, const Scalar &scalar) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        xsimd::batch<T> scalarBatch = xsimd::broadcast(scalar);
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::store_aligned(&result[i], lBatch - scalarBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] - scalar;
        }
        return result;
    }

    /**
    * Applies the Operation - to a Container and a Scalar.
    * @example {1, 2, 3} * 2 = {2, 4, 6}
    * @tparam Container
    * @tparam Scalar
    * @param lhs - multiplicand
    * @param scalar - multiplicand
    * @return a Container
    */
    template<typename Container, typename Scalar>
    inline Container operator*(const Container &lhs, const Scalar &scalar) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        xsimd::batch<T> scalarBatch = xsimd::broadcast(scalar);
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::store_aligned(&result[i], lBatch * scalarBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] * scalar;
        }
        return result;
    }

    /**
     * Applies the Operation / to a Container and a Scalar.
     * @example {2, 4, 6} / 2 = {1, 2, 3}
     * @tparam Container
     * @tparam Scalar
     * @param lhs - the dividend
     * @param scalar - the divisor
     * @return a Container
     */
    template<typename Container, typename Scalar>
    inline Container operator/(const Container &lhs, const Scalar &scalar) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(lhs) - (std::size(lhs) % simd_size);
        Container result{lhs};
        xsimd::batch<T> scalarBatch = xsimd::broadcast(scalar);
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> lBatch = xsimd::load_unaligned(&lhs[i]);
            xsimd::store_aligned(&result[i], lBatch / scalarBatch);
        }
        for (size_t i = ub; i < std::size(lhs); ++i) {
            result[i] = lhs[i] / scalar;
        }
        return result;
    }

    /**
     * Applies the euclidean norm/ L2-norm to a Container (e.g. a vector)
     * @tparam Container - must be iterable
     * @param container - e.g. a vector
     * @return an double containing the L2 norm
     */
    template<typename Container>
    inline double euclideanNorm(const Container &container) {
        typedef typename Container::value_type T;
        constexpr std::size_t simd_size = xsimd::simd_type<T>::size;
        const int ub = std::size(container) - (std::size(container) % simd_size);
        T result{0};
        for (size_t i = 0; i < ub; i+= simd_size) {
            xsimd::batch<T> reg = xsimd::load_unaligned(&container[i]);
            result += xsimd::hadd(reg * reg);
        }
        for (size_t i = ub; i < std::size(container); ++i) {
            result += container[i] * container[i];
        }
        return std::sqrt(result);
    }

    /**
     * Computes the absolute value for each value in the given container
     * @tparam Container - a iterable container, containing numerical values
     * @param container - the container
     * @return a container with the modified values
     */
    template<typename Container>
    inline Container abs(const Container &container) {
        Container ret = container;
        std::transform(std::begin(container), std::end(container), std::begin(ret),
                       [](const auto &element) { return std::abs(element); });
        return ret;
    }

    /**
     * Computes the determinant with the Sarrus rule for a 3x3 matrix.
     * Notice that for square matrices det(A) = det(A^T).
     * @tparam T - a numerical value
     * @param matrix - the 3x3 matrix
     * @return the determinant
     */
    template<typename T>
    inline T det(const Matrix<T, 3, 3> &matrix) {
        return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0]
               + matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0]
               - matrix[0][0] * matrix[1][2] * matrix[2][1] - matrix[0][1] * matrix[1][0] * matrix[2][2];
    }

    /**
     * Computes the transposed of a mxn matrix.
     * @tparam T - the type of the matrix elements
     * @tparam M - the row number
     * @tparam N - the column number
     * @param matrix - the matrix to transpose
     * @return the transposed
     */
    template<typename T, size_t M, size_t N>
    inline Matrix<T, M, N> transpose(const Matrix<T, M, N> &matrix) {
        Matrix<T, N, M> transposed;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                transposed[i][j] = matrix[j][i];
            }
        }
        return transposed;
    }

    /**
    * Returns the cross product of two cartesian vectors.
    * @tparam T - a number
    * @param lhs - left vector
    * @param rhs - right vector
    * @return cross product
    */
    template<typename T>
    inline std::array<T, 3> cross(const std::array<T, 3> &lhs, const std::array<T, 3> &rhs) {
        std::array<T, 3> result{};
        result[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
        result[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
        result[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
        return result;
    }

    /**
    * Returns the dot product of two cartesian vectors.
    * @tparam T - a number
    * @param lhs - left vector
    * @param rhs - right vector
    * @return dot product
    */
    template<typename T>
    inline T dot(const std::array<T, 3> &lhs, const std::array<T, 3> &rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    }

    /**
     * Implements the signum function with a certain EPSILON to absorb rounding errors.
     * @tparam T - a numerical (floating point) value
     * @param val - the value itself
     * @param cutoffEpsilon - the cut-off radius around zero to return 0
     * @return -1, 0, 1 depending on the sign an the given EPSILON
     */
    template<typename T>
    inline int sgn(T val, double cutoffEpsilon) {
        return val < -cutoffEpsilon ? -1 : val > cutoffEpsilon ? 1 : 0;
    }

    /**
     * Concatenates two std::array of different sizes to one array.
     * @tparam T - the shared type of the arrays
     * @tparam M - the size of the first container
     * @tparam N  - the size of the second container
     * @param first - the first array
     * @param second - the second array
     * @return a new array of size M+N with type T
     */
    template<typename T, size_t M, size_t N>
    inline std::array<T, M+N> concat(const std::array<T, M> &first, const std::array<T, N> &second) {
        std::array<T, M+N> result{};
        size_t index = 0;
        for (const auto &el : first) {
            result[index++] = el;
        }
        for (const auto &el : second) {
            result[index++] = el;
        }
        return result;
    }

    /**
     * Operator << for an array of any size.
     * @tparam T - type of the array, must have an << operator overload
     * @tparam N - size of the array
     * @param os - the ostream
     * @param array - the array itself
     * @return ostream
     */
    template<typename T, size_t N>
    inline std::ostream &operator<<(std::ostream &os, const std::array<T, N> &array) {
        os << "[";
        auto it = array.begin();
        auto end = array.end() - 1;
        while (it != end) {
            os << *it << " ";
            ++it;
        }
        os << *it << "]";
        return os;
    }

    template<typename T>
    struct is_stdarray : std::false_type {
    };

    template<typename T, std::size_t N>
    struct is_stdarray<std::array<T, N>> : std::true_type {
    };

}