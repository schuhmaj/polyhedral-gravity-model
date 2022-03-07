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
                    {4, 2, 1},
                    {4, 3, 2},
                    {6, 5, 1},
                    {2, 6, 1},
                    {5, 8, 1},
                    {8, 4, 1},
                    {7, 6, 2},
                    {3, 7, 2},
                    {8, 7, 4},
                    {4, 7, 3},
                    {6, 7, 5},
                    {7, 8, 5}
            }
    };

    polyhedralGravity::Gravity systemUnderTest{_polyhedron};

};

TEST_F(GravityTest, gij_vectors) {
    std::vector<std::array<std::array<double, 3>, 3>> expectedGij;
    auto actualGij = systemUnderTest.calculateGij();

}

