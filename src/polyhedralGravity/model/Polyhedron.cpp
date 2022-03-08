#include "Polyhedron.h"

namespace polyhedralGravity {

    std::pair<const std::array<double, 3> &, const std::array<double, 3> &>
    Polyhedron::getPolyhedralSegment(size_t p, size_t q) const {
        const auto &plane = _faces.at(p);
        const std::array<double, 3> &node1 = _nodes.at(plane.at(q % 3));
        const std::array<double, 3> &node2 = _nodes.at(plane.at((q + 1) % 3));
        return {node1, node2};
    }
}
