#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/calculation/Gravity.h"
#include "polyhedralGravity/model/Polyhedron.h"

/**
 * Contains Tests based on the example from Tsoulis FORTRAN implementation.
 * Hardcoded values taken from his implementation's results.
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
                    {0,   1, 3},
                    {1, 2, 3},
                    {0, 4,  5},
                    {0,   5,  1},
                    {0,   7, 4},
                    {0, 3, 7},
                    {1, 5,  6},
                    {1,   6,  2},
                    {3, 6, 7},
                    {2, 6, 3},
                    {4, 6, 5},
                    {4, 7, 6}
            }
    };

    polyhedralGravity::Gravity systemUnderTest{_polyhedron};

    std::vector<std::array<std::array<double, 3>, 3>> expectedGij{
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, 0.0}, {-20.0, 10.0, 0.0}, {0.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {-20.0, 0.0, 0.0}, {20.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {20.0, 0.0, 0.0}, {-20.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, -10.0}, {0.0, 0.0, 10.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, -10.0}, {0.0, -10.0, 0.0}, {0.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {0.0, 0.0, -10.0}, {0.0, -10.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {0.0, 10.0, 0.0}, {0.0, -10.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, -10.0}, {0.0, 0.0, 10.0}, {0.0, -10.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 0.0, -10.0}, {-20.0, 0.0, 0.0}, {0.0, 0.0, 10.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, -10.0}, {-20.0, 0.0, 10.0}, {20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{20.0, 10.0, 0.0}, {0.0, -10.0, 0.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {20.0, 0.0, 0.0}, {-20.0, -10.0, 0.0}}}
    };

    std::vector<std::array<double, 3>> expectedPlaneUnitNormals{
            {0.0,  -0.0, 1.0},
            {0.0,  -0.0, 1.0},
            {0.0,  -1.0, 0.0},
            {0.0,  -1.0, 0.0},
            {-1.0, -0.0, -0.0},
            {-1.0, 0.0,  0.0},
            {1.0,  -0.0, 0.0},
            {1.0,  -0.0, 0.0},
            {0.0,  1.0,  0.0},
            {0.0,  1.0,  0.0},
            {0.0,  0.0,  -1.0},
            {0.0,  0.0,  -1.0}
    };

    std::vector<std::array<std::array<double, 3>, 3>> expectedSegmentUnitNormals{
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -1.0, -0.0}, {0.4472135954999579, 0.8944271909999159, 0.0}, {-1.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{1.0, 0.0, -0.0}, {0.0, 1.0, 0.0}, {-0.4472135954999579, -0.8944271909999159, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-1.0, -0.0, -0.0}, {0.0, 0.0, -1.0}, {0.4472135954999579, 0.0, 0.8944271909999159}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.4472135954999579, -0.0, -0.8944271909999159}, {1.0, 0.0, -0.0}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.0, 0.7071067811865475, 0.7071067811865475}, {0.0, 0.0, -1.0}, {0.0, -1.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -0.0, 1.0}, {0.0, 1.0, 0.0}, {-0.0, -0.7071067811865475, -0.7071067811865475}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -1.0, -0.0}, {0.0, 0.0, -1.0}, {0.0, 0.7071067811865475, 0.7071067811865475}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.0, -0.7071067811865475, -0.7071067811865475}, {0.0, 1.0, -0.0}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{0.4472135954999579, -0.0, 0.8944271909999159}, {0.0, 0.0, -1.0}, {-1.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{1.0, -0.0, 0.0}, {-0.4472135954999579, 0.0, -0.8944271909999159}, {0.0, 0.0, 1.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-0.4472135954999579, 0.8944271909999159, 0.0}, {1.0, 0.0, 0.0}, {-0.0, -1.0, -0.0}}},
            std::array<std::array<double, 3>, 3>{
                    {{-1.0, 0.0, 0.0}, {-0.0, 1.0, 0.0}, {0.4472135954999579, -0.8944271909999159, 0.0}}}
    };

    std::vector<double> expectedPlaneNormalOrientations{1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0};

    std::vector<polyhedralGravity::HessianPlane> expectedHessianPlanes{
            {0.0,    0.0,    200.0,  -5000.0},
            {0.0,    -0.0,   200.0,  -5000.0},
            {0.0,    -200.0, 0.0,    0.0},
            {0.0,    -200.0, 0.0,    0.0},
            {-100.0, 0.0,    0.0,    -2000.0},
            {-100.0, 0.0,    0.0,    -2000.0},
            {100.0,  0.0,    0.0,    -0.0},
            {100.0,  -0.0,   0.0,    0.0},
            {0.0,    200.0,  0.0,    -2000.0},
            {0.0,    200.0,  0.0,    -2000.0},
            {0.0,    0.0,    -200.0, 3000.0},
            {0.0,    0.0,    -200.0, 3000.0}
    };

    std::vector<double> expectedPlaneDistances{25.0, 25.0, 0.0, 0.0, 20.0, 20.0, 0.0, 0.0, 10.0, 10.0, 15.0, 15.0};

    std::vector<std::array<double, 3>> expectedOrthogonalProjectionPointsOnPlane{
            {0.0,   0.0,  25.0},
            {0.0,   0.0,  25.0},
            {0.0,   0.0,  0.0},
            {0.0,   0.0,  0.0},
            {-20.0, 0.0,  0.0},
            {-20.0, 0.0,  0.0},
            {0.0,   0.0,  0.0},
            {0.0,   0.0,  0.0},
            {0.0,   10.0, 0.0},
            {0.0,   10.0, 0.0},
            {0.0,   0.0,  15.0},
            {0.0,   0.0,  15.0}
    };

    std::vector<std::array<double, 3>> expectedSegmentNormalOrientations{
            {0.0,  0.0,  1.0},
            {0.0,  1.0,  0.0},
            {1.0,  -1.0, 1.0},
            {-1.0, 0.0,  1.0},
            {1.0,  -1.0, 0.0},
            {1.0,  1.0,  -1.0},
            {0.0,  -1.0, 1.0},
            {-1.0, 1.0,  1.0},
            {1.0,  -1.0, 1.0},
            {0.0,  -1.0, 1.0},
            {1.0,  0.0,  0.0},
            {1.0,  1.0,  -1.0}
    };

    std::vector<std::array<std::array<double, 3>, 3>> expectedOrthogonalProjectionPointsOnSegment{
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 25.0}, {0.0, 0.0, 25.0}, {-20.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 25.0}, {-0.0, 10.0, 25.0}, {0.0, 0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, -0.0}, {-0.0, -0.0, 15.0}, {6.0, -0.0, 12.0}}},
            std::array<std::array<double, 3>, 3>{{{6.0, -0.0, 12.0}, {0.0, 0.0, 0.0}, {-0.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, 12.5, 12.5}, {-20.0, -0.0, 15.0}, {-20.0, 0.0, 0.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, 25.0}, {-20.0, 10.0, -0.0}, {-20.0, 12.5, 12.5}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 0.0, 0.0}, {-0.0, -0.0, 15.0}, {-0.0, 12.5, 12.5}}},
            std::array<std::array<double, 3>, 3>{{{-0.0, 12.5, 12.5}, {-0.0, 10.0, -0.0}, {-0.0, -0.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{6.0, 10.0, 12.0}, {-0.0, 10.0, 15.0}, {-20.0, 10.0, -0.0}}},
            std::array<std::array<double, 3>, 3>{{{0.0, 10.0, 0.0}, {6.0, 10.0, 12.0}, {-0.0, 10.0, 25.0}}},
            std::array<std::array<double, 3>, 3>{{{-4.0, 8.0, 15.0}, {0.0, 0.0, 15.0}, {0.0, 0.0, 15.0}}},
            std::array<std::array<double, 3>, 3>{{{-20.0, -0.0, 15.0}, {-0.0, 10.0, 15.0}, {-4.0, 8.0, 15.0}}}
    };

    std::vector<std::array<double, 3>> expectedSegmentDistances{
            {0.0,                0.0,                20.0},
            {0.0,                10.0,               0.0},
            {20.0,               15.0,               13.416407864998739},
            {13.416407864998739, 0.0,                25.0},
            {17.67766952966369,  15.0,               0.0},
            {25.0,               10.0,               17.67766952966369},
            {0.0,                15.0,               17.67766952966369},
            {17.67766952966369,  10.0,               25.0},
            {13.416407864998739, 15.0,               20.0},
            {0.0,                13.416407864998739, 25.0},
            {8.94427190999916,   0.0,                0.0},
            {20.0,               10.0,               8.94427190999916}
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

TEST_F(GravityTest, PlaneNormalOrientations) {
    using namespace testing;

    auto actualPlaneNormalOrientations = systemUnderTest.calculatePlaneNormalOrientations(expectedPlaneUnitNormals);

    ASSERT_THAT(actualPlaneNormalOrientations, ContainerEq(expectedPlaneNormalOrientations));
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

    auto actualHessianPlane = systemUnderTest.calculateFacesToHessianPlanes();

    ASSERT_EQ(actualHessianPlane, expectedHessianPlanes);
}

TEST_F(GravityTest, PlaneDistances) {
    using namespace testing;

    auto actualPlaneDistances = systemUnderTest.calculatePlaneDistances(expectedHessianPlanes);

    ASSERT_THAT(actualPlaneDistances, ContainerEq(expectedPlaneDistances));
}

TEST_F(GravityTest, OrthogonalProjectionPointsOnPlane) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnPlane = systemUnderTest.calculateOrthogonalProjectionPointsOnPlane(
            expectedHessianPlanes, expectedPlaneUnitNormals, expectedPlaneDistances);

    ASSERT_THAT(actualOrthogonalProjectionPointsOnPlane, ContainerEq(expectedOrthogonalProjectionPointsOnPlane));
}

TEST_F(GravityTest, SegmentNormalOrientations) {
    using namespace testing;

    auto actualSegmentNormalOrientations =
            systemUnderTest.calculateSegmentNormalOrientations(expectedSegmentUnitNormals,
                                                               expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_THAT(actualSegmentNormalOrientations, ContainerEq(expectedSegmentNormalOrientations));
}

TEST_F(GravityTest, OrthogonalProjectionPointsOnSegment) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnSegment =
            systemUnderTest.calculateOrthogonalProjectionPointsOnSegments(expectedOrthogonalProjectionPointsOnPlane,
                                                                          expectedSegmentNormalOrientations);

    ASSERT_THAT(actualOrthogonalProjectionPointsOnSegment, ContainerEq(expectedOrthogonalProjectionPointsOnSegment));
}

TEST_F(GravityTest, SegmentDistances) {
    using namespace testing;

    auto actualSegmentDistances =
            systemUnderTest.calculateSegmentDistances(expectedOrthogonalProjectionPointsOnPlane,
                                                      expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualSegmentDistances, ContainerEq(expectedSegmentDistances));
}
