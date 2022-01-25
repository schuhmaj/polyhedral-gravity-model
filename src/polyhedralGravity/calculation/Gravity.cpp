#include "Gravity.h"

void Gravity::calculate() {
    SPDLOG_INFO("Calculate...");
    auto x = calculateGij();
}

std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
    using util::operator-;
    std::vector<std::array<std::array<double, 3>, 3>> g;
    g.reserve(_polyhedron.size());
    for (auto &plane : _polyhedron.getPolyhedralSegments()) {
        std::array<std::array<double, 3>, 3> gi{};
        for (int j = 0; j < gi.size(); ++j) {
            gi[j] = plane.at(j).second - plane.at(j).first;
        }
        g.push_back(gi);
    }
    return g;
}


