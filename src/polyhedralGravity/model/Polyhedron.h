#pragma once

#include <utility>
#include <vector>

/**
 * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
 */
class Polyhedron {

    /**
     * A vector containing the nodes of the polyhedron.
     * Each node is an array of size three containing the xyz coordinates.
     * TODO Explicit numbering scheme?
     */
    std::vector<std::array<double, 3>> nodes;

    /**
     * A vector containing the faces (triangles) of the polyhedron.
     * Each face is an array of size three containing the indices of the nodes forming the face.
     * TODO Explicit numbering scheme?
     */
    std::vector<std::array<size_t, 3>> faces;

    //TODO elements (tetrahedons) required?


public:

    /**
     * Generates a polyhedron from nodes and faces.
     * @param nodes - vector containing the nodes
     * @param faces - vector containing the triangle faces.
     */
    Polyhedron(std::vector<std::array<double, 3>> nodes, std::vector<std::array<size_t, 3>> faces)
            : nodes{std::move(nodes)},
              faces{std::move(faces)} {}

    ~Polyhedron() = default;

};

