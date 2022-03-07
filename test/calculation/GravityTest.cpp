#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/calculation/Gravity.h"
#include "polyhedralGravity/model/Polyhedron.h"

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
                    {3, 1, 0},
                    {3, 2, 1},
                    {5, 4, 0},
                    {1, 5, 0},
                    {4, 7, 0},
                    {7, 3, 0},
                    {6, 5, 1},
                    {2, 6, 1},
                    {7, 6, 3},
                    {3, 6, 2},
                    {5, 6, 4},
                    {6, 7, 4}
            }
    };

    polyhedralGravity::Gravity systemUnderTest{_polyhedron};

};

TEST_F(GravityTest, gij_vectors) {
    using namespace testing;
    using namespace ::polyhedralGravity;

    // G(i, j) = expectedGij(i * 3 + j)
    std::vector<std::array<double, 3>> expectedGij {
            { 20.000000000000000, -10.000000000000000, 0.0000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 10.000000000000000, 0.0000000000000000 },
            { 20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 0.0000000000000000 },
            { -20.000000000000000, 10.000000000000000, 0.0000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, 10.000000000000000 },
            { 20.000000000000000, 0.0000000000000000, -10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, -10.000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 10.000000000000000 },
            { 20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 10.000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, -10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, 10.000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 10.000000000000000, -10.000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, 10.000000000000000 },
            { 0.0000000000000000, 10.000000000000000, -10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, -10.000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 10.000000000000000 },
            { 0.0000000000000000, 10.000000000000000, 0.0000000000000000 },
            { 20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, -10.000000000000000 },
            { 20.000000000000000, 0.0000000000000000, -10.000000000000000 },
            { 0.0000000000000000, 0.0000000000000000, 10.000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, 10.000000000000000, 0.0000000000000000 },
            { -20.000000000000000, -10.000000000000000, 0.0000000000000000 },
            { 20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { -20.000000000000000, 0.0000000000000000, 0.0000000000000000 },
            { 0.0000000000000000, -10.000000000000000, 0.0000000000000000 },
            { 20.000000000000000, 10.000000000000000, 0.0000000000000000 }
    };
    auto actualGij = systemUnderTest.calculateGij();

    for (int i = 0; i < 12; ++i) {
        for(int j = 0; j < 3; ++j) {
            ASSERT_EQ(actualGij.at(i).at(j), expectedGij.at(i * 3 + j)) << "The vector G_ij was not not equal at (i,j)=(" << i
            << ", " << j << ")";
        }
    }
}
