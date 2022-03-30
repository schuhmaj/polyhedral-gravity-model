#pragma once

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <exception>
#include <stdexcept>

namespace polyhedralGravity {

    /**
     * Data structure containing the model data of one polyhedron. This includes nodes, edges (faces) and elements.
     * The index always starts with zero!
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
         *
         * ASSERTS PRE-CONDITION
         * @throws runtime_error if no face contains the node zero indicating mathematical index
         */
        Polyhedron(std::vector<std::array<double, 3>> nodes, std::vector<std::array<size_t, 3>> faces)
                : _nodes{std::move(nodes)},
                  _faces{std::move(faces)} {
            //Checks that the node with index zero is actually used
            if (_faces.end() == std::find_if(_faces.begin(), _faces.end(), [&](auto &face) {
                return face[0] == 0 || face[1] == 0 || face[2] == 0;
            })) {
                throw std::runtime_error("The node with index zero (0) was never used in any face! This is "
                                         "no valid polyhedron. Probable issue: Started counting at one (1).");
            }
        }

        ~Polyhedron() = default;

        [[nodiscard]] const std::vector<std::array<double, 3>> &getNodes() const {
            return _nodes;
        }

        [[nodiscard]] const std::array<double, 3> &getNode(size_t index) const {
            return _nodes[index];
        }

        /**
         * The number of points (nodes) that make up the polyhedron.
         * @return a size_t
         */
        [[nodiscard]] size_t countNodes() const {
            return _nodes.size();
        }

        /**
         * Returns the number of faces (triangles) that make up the polyhedral.
         * @return a size_t
         */
        [[nodiscard]] size_t countFaces() const {
            return _faces.size();
        }

        [[nodiscard]] const std::vector<std::array<size_t, 3>> &getFaces() const {
            return _faces;
        }

    };

}
