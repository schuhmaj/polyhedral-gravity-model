#include "Gravity.h"

namespace polyhedralGravity {

    void Gravity::calculate() {
        SPDLOG_INFO("Calculate...");
        auto x = calculateGij();
    }

    std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
        using util::operator-;
        std::vector<std::array<std::array<double, 3>, 3>> g{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto& face) -> std::array<std::array<double, 3>, 3> {
            const auto &node0 = _polyhedron.getNode(face[0]);
            const auto &node1 = _polyhedron.getNode(face[1]);
            const auto &node2 = _polyhedron.getNode(face[2]);
            return {node1 - node0, node2 - node1, node0 - node2};
        });
        return g;
    }

    std::vector<std::array<double, 3>> Gravity::calculatePlaneUnitNormals() {
        std::vector<std::array<double, 3>> N{_polyhedron.countFaces()};

        return N;
    }

    std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateSegmentUnitNormals() {
        std::vector<std::array<std::array<double, 3>, 3>> n{_polyhedron.countFaces()};

        return n;
    }

}
