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
     * Since every face consist of three nodes, every face consist of three segments. Each segment consists of
     * two nodes.
     * @example face consisting of {1, 2, 3} --> segments: {1, 2}, {2, 3}, {3, 1}
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

    /**
     * The number of points (nodes) that make up the polyhedron.
     * @return a size_t
     */
    [[nodiscard]] size_t size() const {
        return _nodes.size();
    }

    [[nodiscard]] const std::vector<std::array<size_t, 3>> &getFaces() const {
        return _faces;
    }

    /**
     * Returns the endpoints (nodes) of a given polyhedral segment.
     * @param p - the polyhedral face
     * @param q - the index of the segment inside the polyhedral face
     * @return a pair of two 3-dimensional coordinate points
     */
    [[nodiscard]] std::pair<const std::array<double, 3> &, const std::array<double, 3> &>
    getPolyhedralSegment(size_t p, size_t q) const;


    /**
     * Returns a vector consisting of the three segments making up a plane/ face.
     * Every segment is a pair consisting of references to the two cartesian points defining it.
     * @return vector of segments per plane
     */
    [[nodiscard]] std::vector<std::array<std::pair<const std::array<double, 3> &, const std::array<double, 3> &>, 3>>
    getPolyhedralSegments() const;

};

