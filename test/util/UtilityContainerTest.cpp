#include "gtest/gtest.h"

#include <utility>
#include <array>
#include "polyhedralGravity/util/UtilityContainer.h"


TEST(UtilityContainer, VectorContainerPlus) {
    using namespace ::util;
    std::array<int, 3> a{3, 4, 5};
    std::array<int, 3> b{6, 8, 10};

    std::array<int, 3> expectedResult{9, 12, 15};
    auto actualResult = a + b;
    ASSERT_EQ(actualResult, expectedResult);
}

TEST(UtilityContainer, VectorScalarPlus) {
    using namespace ::util;
    std::array<int, 3> a{3, 4, 5};
    int b = 100;

    std::array<int, 3> expectedResult{103, 104, 105};
    auto actualResult = a + b;
    ASSERT_EQ(actualResult, expectedResult);
}