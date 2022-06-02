#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <array>
#include <utility>
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"


/**
 * Contains Tests how the calculation handles a cubic polyhedron
 */
class GravityModelCubeTest : public ::testing::Test {

protected:

    static constexpr double LOCAL_TEST_EPSILON = 10e-10;


    const polyhedralGravity::Polyhedron _cube{{
                                                {-1.0, -1.0, -1.0},
                                                {1.0, -1.0, -1.0},
                                                {1.0, 1.0, -1.0},
                                                {-1.0, 1.0, -1.0},
                                                {-1.0, -1.0, 1.0},
                                                {1.0, -1.0, 1.0},
                                                {1.0, 1.0, 1.0},
                                                {-1.0, 1.0, 1.0}},
                                        {
                                                {1,    3,    2},
                                                {0,   3,    1},
                                                {0,   1,   5},
                                                {0,    5,   4},
                                                {0,    7,    3},
                                                {0,   4,    7},
                                                {1,   2,   6},
                                                {1,    6,   5},
                                                {2, 3, 6},
                                                {3, 7, 6},
                                                {4, 5, 6},
                                                {4, 6, 7}}
    };

    const double _density = 2670.0;

};

TEST_F(GravityModelCubeTest, Origin) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {0.0, 0.0, 0.0};

    double expectedPotential = 1.6965555143811e-6;
    std::array<double, 3> expectedAcceleration = {0.0, 0.0, 0.0};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleEq(), expectedAcceleration));
}

TEST_F(GravityModelCubeTest, Corner01) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {1.0, 1.0, 1.0};

    double expectedPotential = 8.482768661715e-7;
    std::array<double, 3> expectedAcceleration = {3.4540879555272017e-7, 3.4540879555272017e-7, 3.4540879555272017e-7};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}

TEST_F(GravityModelCubeTest, Corner02) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {-1.0, -1.0, -1.0};

    double expectedPotential = 8.482768661715e-7;
    std::array<double, 3> expectedAcceleration = {-3.4540879555272017e-7, -3.4540879555272017e-7, -3.4540879555272017e-7};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}

TEST_F(GravityModelCubeTest, Plane01) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {1.0, 0.0, 0.0};

    double expectedPotential = 1.2779422904244e-6;
    std::array<double, 3> expectedAcceleration = {9.2531666422050232e-7, 0.0, 0.0};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}

TEST_F(GravityModelCubeTest, Plane02) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {0.0, 1.0, 0.0};

    double expectedPotential = 1.2779422904244e-6;
    std::array<double, 3> expectedAcceleration = {0.0, 9.2531666422050232e-7, 0.0};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}

TEST_F(GravityModelCubeTest, Plane03) {
    using namespace testing;
    using namespace polyhedralGravity;

    std::array<double, 3> computationPoint = {0.0, 0.0, 1.0};

    double expectedPotential = 1.2779422904244e-6;
    std::array<double, 3> expectedAcceleration = {0.0, 0.0, 9.2531666422050232e-7};


    GravityModelResult actualResult = GravityModel::evaluate(_cube, _density, computationPoint);

    ASSERT_NEAR(actualResult.gravitationalPotential, expectedPotential, LOCAL_TEST_EPSILON);

    ASSERT_THAT(actualResult.acceleration,
                Pointwise(DoubleNear(LOCAL_TEST_EPSILON), expectedAcceleration));
}