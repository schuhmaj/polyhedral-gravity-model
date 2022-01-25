#pragma once

#include <utility>
#include <vector>
#include <array>

/**
 * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
 */
class Polyhedron {

    /**
     * A vector containing the nodes of the polyhedron.
     * Each node is an array of size three containing the xyz coordinates.
     */
    const std::vector<std::array<double, 3>> _nodes;

    /**
     * A vector containing the faces (triangles) of the polyhedron.
     * Each face is an array of size three containing the indices of the nodes forming the face.
     */
    const std::vector<std::array<size_t, 3>> _faces;


public:

    /**
     * Generates an empty polyhedron.
     */
    Polyhedron()
            : _nodes{},
              _faces{} {}

    /**
     * Generates a polyhedron from nodes and faces.
     * @param nodes - vector containing the nodes
     * @param faces - vector containing the triangle faces.
     */
    Polyhedron(std::vector<std::array<double, 3>> nodes, std::vector<std::array<size_t, 3>> faces)
            : _nodes{std::move(nodes)},
              _faces{std::move(faces)} {}

    ~Polyhedron() = default;

    [[nodiscard]] const std::vector<std::array<double, 3>> &getNodes() const {
        return _nodes;
    }

    [[nodiscard]] const std::vector<std::array<size_t, 3>> &getFaces() const {
        return _faces;
    }

};

