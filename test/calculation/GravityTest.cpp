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
    std::vector<std::array<std::array<double, 3>, 3>> expectedGij;
    auto actualGij = systemUnderTest.calculateGij();

}

