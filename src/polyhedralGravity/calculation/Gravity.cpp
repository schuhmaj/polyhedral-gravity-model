#include "Gravity.h"

namespace polyhedralGravity {

    void Gravity::calculate() {
        SPDLOG_INFO("Calculate...");
        auto x = calculateGij();
    }

    std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
        using util::operator-;
        std::vector<std::array<std::array<double, 3>, 3>> g;
        g.resize(_polyhedron.countFaces() * 3);
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto& face) -> std::array<std::array<double, 3>, 3> {
            const auto &node0 = _polyhedron.getNode(face[0]);
            const auto &node1 = _polyhedron.getNode(face[1]);
            const auto &node2 = _polyhedron.getNode(face[2]);
            return {node1 - node0, node2 - node1, node0 - node2};
        });
        return g;
    }

}
