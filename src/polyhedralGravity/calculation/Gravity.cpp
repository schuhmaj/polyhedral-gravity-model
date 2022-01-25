#include "Gravity.h"

void Gravity::calculate() {
    SPDLOG_INFO("Calculate...");
    auto x = calculateGij();
}

std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
    using util::operator-;
    std::vector<std::array<std::array<double, 3>, 3>> g;
    g.reserve(_polyhedron.size());
    for (auto &face: _polyhedron.getFaces()) {
        auto &node0 = _polyhedron.getNode(face[0]);
        auto &node1 = _polyhedron.getNode(face[1]);
        auto &node2 = _polyhedron.getNode(face[2]);
        g.push_back({
                            node1 - node0,
                            node2 - node1,
                            node0 - node2
                    });
    }
    return g;
}


