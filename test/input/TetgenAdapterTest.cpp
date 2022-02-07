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
    };

    std::vector<std::string> simpleFiles;
    simpleFiles.emplace_back("resources/TetgenAdapterTestReadSimple.node");
    simpleFiles.emplace_back("resources/TetgenAdapterTestReadSimple.face");

    TetgenAdapter tetgenAdapter{simpleFiles};
    auto actualPolyhedron = tetgenAdapter.getPolyhedron();

    ASSERT_THAT(actualPolyhedron.getFaces(), ContainerEq(expectedFaces));
}