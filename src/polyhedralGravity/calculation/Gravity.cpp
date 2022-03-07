#include "Gravity.h"

namespace polyhedralGravity {

    void Gravity::calculate() {
        SPDLOG_INFO("Calculate...");
        auto x = calculateGij();
    }

    std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
        std::vector<std::array<std::array<double, 3>, 3>> g{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto &face) -> std::array<std::array<double, 3>, 3> {
                           using util::operator-;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           return {node1 - node0, node2 - node1, node0 - node2};
                       });
        return g;
    }

    std::vector<std::array<double, 3>>
    Gravity::calculatePlaneUnitNormals(const std::vector<std::array<std::array<double, 3>, 3>> &g) {
        std::vector<std::array<double, 3>> planeUnitNormal{_polyhedron.countFaces()};
        std::transform(g.cbegin(), g.cend(), planeUnitNormal.begin(), [](const auto &gi) -> std::array<double, 3> {
            using namespace util;
            const std::array<double, 3> crossProduct = cross(gi[0], gi[1]);
            const double norm = euclideanNorm(crossProduct);
            return crossProduct / norm;
        });
        return planeUnitNormal;
    }

    std::vector<std::array<std::array<double, 3>, 3>>
    Gravity::calculateSegmentUnitNormals(const std::vector<std::array<std::array<double, 3>, 3>> &g,
                                         const std::vector<std::array<double, 3>> &planeUnitNormals) {
        std::vector<std::array<std::array<double, 3>, 3>> segmentUnitNormal{_polyhedron.countFaces()};

        return segmentUnitNormal;
    }

}
