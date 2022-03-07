#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include "polyhedralGravity/input/TetgenAdapter.h"

class TetgenAdapterTest : public ::testing::Test {

};

TEST_F(TetgenAdapterTest, readSimpleNode) {
    using namespace testing;
    using namespace ::polyhedralGravity;
    std::vector<std::array<double, 3>> expectedNodes = {
            {-20, 0,  25},
            {0,   0,  25},
            {0,   10, 25},
            {-20, 10, 25},
            {-20, 0,  15},
            {0,   0,  15},
            {0,   10, 15},
            {-20, 10, 15}
    };

    std::vector<std::string> simpleFiles {
        "resources/TetgenAdapterTestReadSimple.node",
        "resources/TetgenAdapterTestReadSimple.face",
    };

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    ASSERT_THAT(actualPolyhedron.getNodes(), ContainerEq(expectedNodes));
}


TEST_F(TetgenAdapterTest, readSimpleFace) {
    using namespace testing;
    using namespace ::polyhedralGravity;
    std::vector<std::array<size_t, 3>> expectedFaces = {
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
    };

    std::vector<std::string> simpleFiles;
    simpleFiles.emplace_back("resources/TetgenAdapterTestReadSimple.node");
    simpleFiles.emplace_back("resources/TetgenAdapterTestReadSimple.face");

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    ASSERT_THAT(actualPolyhedron.getFaces(), ContainerEq(expectedFaces));
}