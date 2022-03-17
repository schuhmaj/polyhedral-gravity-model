#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/calculation/Gravity.h"
#include "polyhedralGravity/model/Polyhedron.h"

/**
 * Contains Tests based on the example from Tsoulis FORTRAN implementation.
 * Harcoded values taken from his implementation's results.
 */
class GravityTest : public ::testing::Test {

protected:
    //New polyhedron with given vertices and faces
    //this is the base example from Tsoulis
    polyhedralGravity::Polyhedron _polyhedron{
            {
                    {-20, 0, 25},
                    {0, 0, 25},
                    {0, 10, 25},
                    {-20, 10, 25},
                    {-20, 0, 15},
                    {0, 0, 15},
                    {0, 10, 15},
                    {-20, 10, 15}
            },
            {
                    {3,   1, 0},
                    {3, 2, 1},
                    {5, 4,  0},
                    {1,   5,  0},
                    {4,   7, 0},
                    {7, 3, 0},
                    {6, 5,  1},
                    {2,   6,  1},
                    {7, 6, 3},
                    {3, 6, 2},
                    {5, 6, 4},
                    {6, 7, 4}
            }
    };

    polyhedralGravity::Gravity systemUnderTest{_polyhedron};

    std::vector<std::array<std::array<double, 3>, 3>> expectedGij{
            std::array<std::array<double, 3>, 3>{{{20.000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {-20.000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 10.000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{20.000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {-20.000000000000000, 10.000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-20.000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {20.000000000000000, 0.0000000000000000, -10.000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 0.0000000000000000, -10.000000000000000},
                                                  {-20.000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {20.000000000000000, 0.0000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 10.000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, -10.000000000000000, 10.000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, -10.000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {0.0000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 10.000000000000000, -10.000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {0.0000000000000000, 10.000000000000000, -10.000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 0.0000000000000000, -10.000000000000000},
                                                  {0.0000000000000000, -10.000000000000000, 10.000000000000000},
                                                  {0.0000000000000000, 10.000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{20.000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {-20.000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, -10.000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{20.000000000000000, 0.0000000000000000, -10.000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, 10.000000000000000},
                                                  {-20.000000000000000, 0.0000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 10.000000000000000, 0.0000000000000000},
                                                  {-20.000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {20.000000000000000, 0.0000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-20.000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, -10.000000000000000, 0.0000000000000000},
                                                  {20.000000000000000, 10.000000000000000, 0.0000000000000000}}},

    };

    std::vector<std::array<double, 3>> expectedPlaneUnitNormals{
            {-0.0000000000000000, -0.0000000000000000, -1.0000000000000000},
            {0.0000000000000000,  0.0000000000000000,  -1.0000000000000000},
            {0.0000000000000000,  1.0000000000000000,  -0.0000000000000000},
            {0.0000000000000000,  1.0000000000000000,  0.0000000000000000},
            {1.0000000000000000,  0.0000000000000000,  -0.0000000000000000},
            {1.0000000000000000,  0.0000000000000000,  -0.0000000000000000},
            {-1.0000000000000000, 0.0000000000000000,  0.0000000000000000},
            {-1.0000000000000000, -0.0000000000000000, -0.0000000000000000},
            {0.0000000000000000,  -1.0000000000000000, 0.0000000000000000},
            {0.0000000000000000,  -1.0000000000000000, 0.0000000000000000},
            {0.0000000000000000,  -0.0000000000000000, 1.0000000000000000},
            {0.0000000000000000,  0.0000000000000000,  1.0000000000000000}
    };

    std::vector<std::array<std::array<double, 3>, 3>> expectedSegmentUnitNormals{
            std::array<std::array<double, 3>, 3>{{{0.44721359549995793, 0.89442719099991586, -0.0000000000000000},
                                                  {0.0000000000000000, -1.0000000000000000, 0.0000000000000000},
                                                  {-1.0000000000000000, 0.0000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, 1.0000000000000000, 0.0000000000000000},
                                                  {1.0000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {-0.44721359549995793, -0.89442719099991586, -0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, 0.0000000000000000, -1.0000000000000000},
                                                  {-1.0000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.44721359549995793, 0.0000000000000000, 0.89442719099991586}}},
            std::array<std::array<double, 3>, 3>{{{1.0000000000000000, -0.0000000000000000, 0.0000000000000000},
                                                  {-0.44721359549995793, 0.0000000000000000, -0.89442719099991586},
                                                  {0.0000000000000000, 0.0000000000000000, 1.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, 0.0000000000000000, -1.0000000000000000},
                                                  {0.0000000000000000, 0.70710678118654746, 0.70710678118654746},
                                                  {0.0000000000000000, -1.0000000000000000, 0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, 1.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, 1.0000000000000000},
                                                  {0.0000000000000000, -0.70710678118654746, -0.70710678118654746}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, -0.0000000000000000, -1.0000000000000000},
                                                  {0.0000000000000000, -1.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, 0.70710678118654746, 0.70710678118654746}}},
            std::array<std::array<double, 3>, 3>{{{-0.0000000000000000, 1.0000000000000000, 0.0000000000000000},
                                                  {0.0000000000000000, -0.70710678118654746, -0.70710678118654746},
                                                  {0.0000000000000000, 0.0000000000000000, 1.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 0.0000000000000000, -1.0000000000000000},
                                                  {0.44721359549995793, 0.0000000000000000, 0.89442719099991586},
                                                  {-1.0000000000000000, -0.0000000000000000, -0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{-0.44721359549995793, -0.0000000000000000, -0.89442719099991586},
                                                  {1.0000000000000000, 0.0000000000000000, -0.0000000000000000},
                                                  {0.0000000000000000, 0.0000000000000000, 1.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{1.0000000000000000, 0.0000000000000000, -0.0000000000000000},
                                                  {-0.44721359549995793, 0.89442719099991586, 0.0000000000000000},
                                                  {0.0000000000000000, -1.0000000000000000, -0.0000000000000000}}},
            std::array<std::array<double, 3>, 3>{{{0.0000000000000000, 1.0000000000000000, -0.0000000000000000},
                                                  {-1.0000000000000000, 0.0000000000000000, 0.0000000000000000},
                                                  {0.44721359549995793, -0.89442719099991586, 0.0000000000000000}}},

    };

    std::vector<double> expectedSigmaP{-1.0, -1.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, -1.0, -1.0, 1.0, 1.0};

    std::vector<polyhedralGravity::HessianPlane> expectedHessianPlane{
            {0.0000000000000000,  0.0000000000000000,  -200.00000000000000, 5000.0000000000000},
            {0.0000000000000000,  0.0000000000000000,  -200.00000000000000, 5000.0000000000000},
            {0.0000000000000000,  200.00000000000000,  0.0000000000000000,  -0.0000000000000000},
            {0.0000000000000000,  200.00000000000000,  0.0000000000000000,  -0.0000000000000000},
            {100.00000000000000,  0.0000000000000000,  0.0000000000000000,  2000.0000000000000},
            {100.00000000000000,  0.0000000000000000,  -0.0000000000000000, 2000.0000000000000},
            {-100.00000000000000, 0.0000000000000000,  0.0000000000000000,  0.0000000000000000},
            {-100.00000000000000, -0.0000000000000000, -0.0000000000000000, 0.0000000000000000},
            {0.0000000000000000,  -200.00000000000000, 0.0000000000000000,  2000.0000000000000},
            {0.0000000000000000,  -200.00000000000000, 0.0000000000000000,  2000.0000000000000},
            {0.0000000000000000,  -0.0000000000000000, 200.00000000000000,  -3000.0000000000000},
            {0.0000000000000000,  0.0000000000000000,  200.00000000000000,  -3000.0000000000000}
    };

    std::vector<double> expectedPlaneDistance{
            25.000000000000000,
            25.000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            20.000000000000000,
            20.000000000000000,
            0.0000000000000000,
            0.0000000000000000,
            10.000000000000000,
            10.000000000000000,
            15.000000000000000,
            15.000000000000000
    };

};

TEST_F(GravityTest, GijVectors) {
    using namespace testing;

    auto actualGij = systemUnderTest.calculateGij();

    ASSERT_THAT(actualGij, ContainerEq(expectedGij));
}

TEST_F(GravityTest, PlaneUnitNormals) {
    using namespace testing;

    auto actualPlaneUnitNormals = systemUnderTest.calculatePlaneUnitNormals(expectedGij);

    ASSERT_THAT(actualPlaneUnitNormals, ContainerEq(expectedPlaneUnitNormals));
}

TEST_F(GravityTest, SegmentUnitNormals) {
    using namespace testing;

    auto actualSegmentUnitNormals = systemUnderTest.calculateSegmentUnitNormals(expectedGij, expectedPlaneUnitNormals);

    ASSERT_THAT(actualSegmentUnitNormals, ContainerEq(expectedSegmentUnitNormals));
}

TEST_F(GravityTest, SigmaP) {
    using namespace testing;

    auto actualSigmaP = systemUnderTest.calculateSigmaP(expectedPlaneUnitNormals);

    ASSERT_THAT(actualSigmaP, ContainerEq(expectedSigmaP));
}

TEST_F(GravityTest, SimpleHessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    HessianPlane expectedHessian{2, -8, 5, -18};

    auto actualHessianPlane = systemUnderTest.computeHessianPlane({1, -2, 0}, {3, 1, 4}, {0, -1, 2});

    ASSERT_DOUBLE_EQ(actualHessianPlane.a, expectedHessian.a);
    ASSERT_DOUBLE_EQ(actualHessianPlane.b, expectedHessian.b);
    ASSERT_DOUBLE_EQ(actualHessianPlane.c, expectedHessian.c);
    ASSERT_DOUBLE_EQ(actualHessianPlane.d, expectedHessian.d);
}

TEST_F(GravityTest, HessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualHessianPlane = systemUnderTest.calculateFaceToHessianPlane();

    ASSERT_EQ(actualHessianPlane, expectedHessianPlane);
}

TEST_F(GravityTest, PlaneDistances) {
    using namespace testing;

    auto actualPlaneDistances = systemUnderTest.calculatePlaneDistance(expectedHessianPlane);

    ASSERT_THAT(actualPlaneDistances, ContainerEq(expectedPlaneDistance));
}

