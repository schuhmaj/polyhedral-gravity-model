#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <string>
#include "polyhedralGravity/input/TetgenAdapter.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/model/Polyhedron.h"


/**
 * Contains Tests based on the Eros mesh taken from
 * https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes (last accessed: 07.04.2022)
 *
 * The values are in the corresponding files in test/resources which are used to check the C++
 * implementation are calculated by the Tsoulis reference implementation in FORTRAN.
 *
 */
class GravityModelBigTest : public ::testing::Test {

protected:

    const size_t countFaces = 14744;
    const size_t countNodesPerFace = 3;

    polyhedralGravity::Polyhedron _polyhedron{
            polyhedralGravity::TetgenAdapter{
                    {"resources/GravityModelBigTest.node", "resources/GravityModelBigTest.face"}}.getPolyhedron()};

    polyhedralGravity::GravityModel systemUnderTest{_polyhedron};

    std::vector<std::array<std::array<double, 3>, 3>> expectedGij;

    std::vector<std::array<double, 3>> expectedPlaneUnitNormals;

    std::vector<std::array<std::array<double, 3>, 3>> expectedSegmentUnitNormals;

    std::vector<double> expectedPlaneNormalOrientations;

    std::vector<polyhedralGravity::HessianPlane> expectedHessianPlanes;

    std::vector<double> expectedPlaneDistances;

    std::vector<std::array<double, 3>> expectedOrthogonalProjectionPointsOnPlane;

    std::vector<std::array<double, 3>> expectedSegmentNormalOrientations;

    std::vector<std::array<std::array<double, 3>, 3>> expectedOrthogonalProjectionPointsOnSegment;

    std::vector<std::array<double, 3>> expectedSegmentDistances;

    std::vector<std::array<polyhedralGravity::Distance, 3>> expectedDistancesPerSegmentEndpoint;

    std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>> expectedTranscendentalExpressions;


    std::vector<std::pair<double, std::array<double, 3>>> expectedSingularityTerms;

    std::vector<double> expectedAlphaSingularityTerms;

    std::vector<std::array<double, 3>> expectedBetaSingularityTerms;

public:

    [[nodiscard]] std::vector<std::array<std::array<double, 3>, 3>>
    readTwoDimensionalCartesian(const std::string &filename) const {
        std::vector<std::array<std::array<double, 3>, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x, y, z;
            if (!(linestream >> x >> y >> z)) {
                break;
            }
            result[i / 3][i % 3] = std::array<double, 3>{x, y, z};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>> readOneDimensionalCartesian(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x, y, z;
            if (!(linestream >> x >> y >> z)) {
                break;
            }
            result[i] = std::array<double, 3>{x, y, z};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>> readTwoDimensionalValue(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x;
            if (!(linestream >> x)) {
                break;
            }
            result[i / 3][i % 3] = x;
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<double> readOneDimensionalValue(const std::string &filename) const {
        std::vector<double> result(countFaces, 0.0);
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double x;
            if (!(linestream >> x)) {
                break;
            }
            result[i] = x;
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<polyhedralGravity::HessianPlane>
    readHessianPlanes(const std::string &filename) const {
        std::vector<polyhedralGravity::HessianPlane> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double a, b, c, d;
            if (!(linestream >> a >> b >> c >> d)) {
                break;
            }
            result[i] = polyhedralGravity::HessianPlane{a, b, c, d};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<polyhedralGravity::Distance, 3>>
    readDistances(const std::string &filename) const {
        std::vector<std::array<polyhedralGravity::Distance, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double l1, l2, s1, s2;
            if (!(linestream >> l1 >> l2 >> s1 >> s2)) {
                break;
            }
            result[i / 3][i % 3] = polyhedralGravity::Distance{l1, l2, s1, s2};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>>
    readTranscendentalExpressions(const std::string &filename) const {
        std::vector<std::array<polyhedralGravity::TranscendentalExpression, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        int i = 0;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            double ln, an;
            if (!(linestream >> ln >> an)) {
                break;
            }
            result[i / 3][i % 3] = polyhedralGravity::TranscendentalExpression{ln, an};
            i += 1;;
        }
        return result;
    }

    [[nodiscard]] std::vector<std::array<double, 3>>
    readBetaSingularities(const std::string &filename) const {
        std::vector<std::array<double, 3>> result{countFaces};
        std::ifstream infile(filename);
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream linestream(line);
            int i, j;
            double sng;
            if (!(linestream >> j >> i >> sng)) {
                break;
            }
            result[i][j] = sng;
        }
        return result;
    }

    GravityModelBigTest() : ::testing::Test() {
        using namespace polyhedralGravity;
        using namespace util;

        expectedGij = readTwoDimensionalCartesian("resources/GravityModelBigTestExpectedGij.txt");
        expectedPlaneUnitNormals =
                readOneDimensionalCartesian("resources/GravityModelBigTestExpectedPlaneUnitNormals.txt");
        expectedSegmentUnitNormals =
                readTwoDimensionalCartesian("resources/GravityModelBigTestExpectedSegmentUnitNormals.txt");
        expectedPlaneNormalOrientations =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedPlaneOrientation.txt");
        expectedHessianPlanes =
                readHessianPlanes("resources/GravityModelBigTestExpectedHessianPlanes.txt");
        expectedPlaneDistances =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedPlaneDistances.txt");
        expectedOrthogonalProjectionPointsOnPlane =
                readOneDimensionalCartesian("resources/GravityModelBigTestExpectedOrthogonalPlaneProjectionPoints.txt");
        expectedSegmentNormalOrientations =
                readTwoDimensionalValue("resources/GravityModelBigTestExpectedSegmentOrientation.txt");
        expectedOrthogonalProjectionPointsOnSegment = readTwoDimensionalCartesian(
                "resources/GravityModelBigTestExpectedOrthogonalSegmentProjectionPoints.txt");
        expectedSegmentDistances =
                readTwoDimensionalValue("resources/GravityModelBigTestExpectedSegmentDistances.txt");
        expectedDistancesPerSegmentEndpoint =
                readDistances("resources/GravityModelBigTestExpectedDistances.txt");
        expectedTranscendentalExpressions =
                readTranscendentalExpressions("resources/GravityModelBigTestExpectedTranscendentalExpressions.txt");
        expectedAlphaSingularityTerms =
                readOneDimensionalValue("resources/GravityModelBigTestExpectedAlphaSingularities.txt");
        expectedBetaSingularityTerms =
                readBetaSingularities("resources/GravityModelBigTestExpectedBetaSingularities.txt");

        expectedSingularityTerms.resize(expectedAlphaSingularityTerms.size());
        for (int i = 0; i < expectedAlphaSingularityTerms.size(); ++i) {
            expectedSingularityTerms[i] =
                    std::make_pair(expectedAlphaSingularityTerms[i], expectedBetaSingularityTerms[i]);
        }
    }

};

TEST_F(GravityModelBigTest, GijVectors) {
    using namespace testing;

    auto actualGij = systemUnderTest.calculateGij();

    ASSERT_THAT(actualGij, ContainerEq(expectedGij));
}

TEST_F(GravityModelBigTest, PlaneUnitNormals) {
    using namespace testing;

    auto actualPlaneUnitNormals = systemUnderTest.calculatePlaneUnitNormals(expectedGij);

    ASSERT_THAT(actualPlaneUnitNormals, ContainerEq(expectedPlaneUnitNormals));
}

TEST_F(GravityModelBigTest, SegmentUnitNormals) {
    using namespace testing;

    auto actualSegmentUnitNormals = systemUnderTest.calculateSegmentUnitNormals(expectedGij, expectedPlaneUnitNormals);

    ASSERT_THAT(actualSegmentUnitNormals, ContainerEq(expectedSegmentUnitNormals));
}

TEST_F(GravityModelBigTest, PlaneNormalOrientations) {
    using namespace testing;

    auto actualPlaneNormalOrientations = systemUnderTest.calculatePlaneNormalOrientations(expectedPlaneUnitNormals);

    ASSERT_THAT(actualPlaneNormalOrientations, ContainerEq(expectedPlaneNormalOrientations));
}

TEST_F(GravityModelBigTest, HessianPlane) {
    using namespace testing;
    using namespace polyhedralGravity;

    auto actualHessianPlane = systemUnderTest.calculateFacesToHessianPlanes();

    ASSERT_EQ(actualHessianPlane, expectedHessianPlanes);
}

TEST_F(GravityModelBigTest, PlaneDistances) {
    using namespace testing;

    auto actualPlaneDistances = systemUnderTest.calculatePlaneDistances(expectedHessianPlanes);

    ASSERT_THAT(actualPlaneDistances, ContainerEq(expectedPlaneDistances));
}

TEST_F(GravityModelBigTest, OrthogonalProjectionPointsOnPlane) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnPlane = systemUnderTest.calculateOrthogonalProjectionPointsOnPlane(
            expectedHessianPlanes, expectedPlaneUnitNormals, expectedPlaneDistances);

    for (size_t i = 0; i < actualOrthogonalProjectionPointsOnPlane.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(
                    actualOrthogonalProjectionPointsOnPlane[i][j],
                    expectedOrthogonalProjectionPointsOnPlane[i][j])
                    << "Difference for P' of plane=" << i << " and coordinate-Nr.=" << j;
        }
    }

    //ASSERT_THAT(actualOrthogonalProjectionPointsOnPlane, ContainerEq(expectedOrthogonalProjectionPointsOnPlane));
}

TEST_F(GravityModelBigTest, SegmentNormalOrientations) {
    using namespace testing;

    auto actualSegmentNormalOrientations =
            systemUnderTest.calculateSegmentNormalOrientations(expectedSegmentUnitNormals,
                                                               expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_THAT(actualSegmentNormalOrientations, ContainerEq(expectedSegmentNormalOrientations));
}

TEST_F(GravityModelBigTest, OrthogonalProjectionPointsOnSegment) {
    using namespace testing;

    auto actualOrthogonalProjectionPointsOnSegment =
            systemUnderTest.calculateOrthogonalProjectionPointsOnSegments(expectedOrthogonalProjectionPointsOnPlane,
                                                                          expectedSegmentNormalOrientations);

    for (size_t i = 0; i < actualOrthogonalProjectionPointsOnSegment.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                EXPECT_NEAR(
                        actualOrthogonalProjectionPointsOnSegment[i][j][k],
                        expectedOrthogonalProjectionPointsOnSegment[i][j][k], 1e16)
                                    << "Difference for P'' of segment=(" << i << ", " << j << ") and coordinate-Nr."
                                    << k;
            }
        }
    }

    //ASSERT_THAT(actualOrthogonalProjectionPointsOnSegment, ContainerEq(expectedOrthogonalProjectionPointsOnSegment));
}

TEST_F(GravityModelBigTest, SegmentDistances) {
    using namespace testing;

    auto actualSegmentDistances =
            systemUnderTest.calculateSegmentDistances(expectedOrthogonalProjectionPointsOnPlane,
                                                      expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualSegmentDistances, ContainerEq(expectedSegmentDistances));
}

TEST_F(GravityModelBigTest, DistancesPerSegmentEndpoint) {
    using namespace testing;

    auto actualDistancesPerSegmentEndpoint =
            systemUnderTest.calculateDistances(expectedGij, expectedOrthogonalProjectionPointsOnSegment);

    ASSERT_THAT(actualDistancesPerSegmentEndpoint, ContainerEq(expectedDistancesPerSegmentEndpoint));
}

TEST_F(GravityModelBigTest, TranscendentalExpressions) {
    using namespace testing;

    auto actualTranscendentalExpressions =
            systemUnderTest.calculateTranscendentalExpressions(expectedDistancesPerSegmentEndpoint,
                                                               expectedPlaneDistances,
                                                               expectedSegmentDistances,
                                                               expectedSegmentNormalOrientations,
                                                               expectedOrthogonalProjectionPointsOnPlane);

    ASSERT_THAT(actualTranscendentalExpressions, ContainerEq(expectedTranscendentalExpressions));
}

TEST_F(GravityModelBigTest, SingularityTerms) {
    using namespace testing;

    auto actualSingularityTerms =
            systemUnderTest.calculateSingularityTerms(expectedGij, expectedSegmentNormalOrientations,
                                                      expectedOrthogonalProjectionPointsOnPlane,
                                                      expectedPlaneDistances, expectedPlaneNormalOrientations,
                                                      expectedPlaneUnitNormals);

    ASSERT_THAT(actualSingularityTerms, ContainerEq(expectedSingularityTerms));
}
