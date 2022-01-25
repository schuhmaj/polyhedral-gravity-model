#include "Polyhedron.h"

std::pair<const std::array<double, 3> &, const std::array<double, 3> &>
Polyhedron::getPolyhedralSegment(size_t p, size_t q) const {
    const auto &plane = _faces.at(p);
    const std::array<double, 3> &node1 = _nodes.at(plane.at(q % 3));
    const std::array<double, 3> &node2 = _nodes.at(plane.at((q + 1) % 3));
    return {node1, node2};
}

std::vector<std::array<std::pair<const std::array<double, 3> &, const std::array<double, 3> &>, 3>>
Polyhedron::getPolyhedralSegments() const {
    std::vector<std::array<std::pair<const std::array<double, 3> &, const std::array<double, 3> &>, 3>> result{};
    result.reserve(_faces.size());
    for (auto &plane: _faces) {
        const std::array<double, 3> &node1 = _nodes.at(plane.at(0));
        const std::array<double, 3> &node2 = _nodes.at(plane.at(1));
        const std::array<double, 3> &node3 = _nodes.at(plane.at(2));
        result.push_back({
                                 std::make_pair(node1, node2),
                                 std::make_pair(node2, node3),
                                 std::make_pair(node3, node1)
                         });
    }
    return result;
}
