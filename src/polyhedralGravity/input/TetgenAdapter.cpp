#include "TetgenAdapter.h"
#include "spdlog/spdlog.h"

namespace polyhedralGravity {

    Polyhedron TetgenAdapter::getPolyhedron() {
        //0. Step: Reset the internal state from previous reads
        _hasNodes = false;
        _hasFaces = false;
        _hasElements = false;

        //1. Step: Read in from files
        for (const auto &fileName: _fileNames) {
            size_t pos = fileName.find_last_of('.');
            std::string name = fileName.substr(0, pos);
            std::string suffix = fileName.substr(pos + 1);
            _suffixToOperation.at(suffix)(name);
        }

        //2. Convert tetgenio to Polyhedron
        return convertTetgenToPolyhedron();
    }

    void TetgenAdapter::readNode(const std::string &filename) {
        if (!_hasNodes) {
            try {
                _tetgenio.load_node(const_cast<char *>(filename.c_str()));
                _hasNodes = true;
            } catch (...) {
                throw std::runtime_error(
                        "The nodes were not read because of an error in Tetgen!"
                        );
            }
        } else {
            throw std::runtime_error(
                    "The Polyhedron already has well defined nodes! The information of " + filename
                    + ".node is redundant!");
        }
    }

    void TetgenAdapter::readFace(const std::string &filename) {
        if (!_hasFaces) {
            try {
            _tetgenio.load_face(const_cast<char *>(filename.c_str()));
            _hasFaces = true;
            } catch (...) {
                throw std::runtime_error(
                        "The faces were read unsuccessfully cause of an fault in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        } else {
            throw std::runtime_error(
                    "The Polyhedron already has well defined faces! The information of " + filename
                    + ".face is redundant!");
        }
    }

    void TetgenAdapter::readElements(const std::string &filename) {
        if (!_hasElements) {
            try {
            _tetgenio.load_elem(const_cast<char *>(filename.c_str()));
            _hasElements = true;
            } catch (...) {
                throw std::runtime_error(
                        "The elements were not read because of an error in Tetgen!"
                );
            }
        } else {
            throw std::runtime_error(
                    "The Polyhedron already has well defined elements! The information of " + filename
                    + ".ele is redundant!");
        }
    }

    Polyhedron TetgenAdapter::convertTetgenToPolyhedron() const {
        std::vector<std::array<double, 3>> nodes{};
        nodes.reserve(_tetgenio.numberofpoints);
        for (size_t i = 0; i < _tetgenio.numberofpoints * 3; i += 3) {
            nodes.push_back({_tetgenio.pointlist[i],
                             _tetgenio.pointlist[i + 1],
                             _tetgenio.pointlist[i + 2]});
        }

        std::vector<std::array<size_t, 3>> faces{};
        faces.reserve(_tetgenio.numberoftrifaces);
        for (size_t i = 0; i < _tetgenio.numberoftrifaces * 3; i += 3) {
            faces.push_back({static_cast<size_t>(_tetgenio.trifacelist[i]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 1]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 2])});
        }
        return {nodes, faces};
    }

}
