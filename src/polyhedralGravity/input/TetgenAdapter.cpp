#include "TetgenAdapter.h"
#include "spdlog/spdlog.h"

namespace polyhedralGravity {

    Polyhedron TetgenAdapter::getPolyhedron() {
        //1. Step: Read in from files
        for (const auto &fileName: _fileNames) {
            size_t pos = fileName.find_last_of('.');
            std::string name = fileName.substr(0, pos);
            std::string suffix = fileName.substr(pos + 1);
            _suffixToOperation.at(suffix)(name);
        }

        //2. Convert tetgenio to Polyhedron
        return {_vertices, _faces};
    }

    void TetgenAdapter::readNode(const std::string &filename) {
        if (this->checkIntegrity(filename, 'v')) {
            try {
                _tetgenio.load_node(const_cast<char *>(filename.c_str()));
                this->addVertices();
            } catch (...) {
                throw std::runtime_error(
                        "The nodes were not read because of an error in Tetgen!"
                );
            }
        }
    }

    void TetgenAdapter::readFace(const std::string &filename) {
        if (this->checkIntegrity(filename, 'f')) {
            try {
                _tetgenio.load_face(const_cast<char *>(filename.c_str()));
                this->addFacesByTrifaces();
            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    void TetgenAdapter::readOff(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_off(const_cast<char *>(filename.c_str()));
                this->addVertices();
                this->addFacesByFacetList();
            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    void TetgenAdapter::readPly(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_ply(const_cast<char *>(filename.c_str()));
                this->addVertices();
                this->addFacesByFacetList();
            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    void TetgenAdapter::readStl(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_stl(const_cast<char *>(filename.c_str()));
                this->addVertices();
                this->addFacesByFacetList();

                tetgenbehavior tetgenbehavior;
                tetgenbehavior.zeroindex = 1;
                tetrahedralize(&tetgenbehavior, &_tetgenio, &_tetgenio);

                this->addVertices();
                this->addFacesByTrifaces();
            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    void TetgenAdapter::readMesh(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_medit(const_cast<char *>(filename.c_str()), 0);
                this->addVertices();
                this->addFacesByFacetList();
                // Additionally .mesh files start counting the index at 1, instead of 0, so decrement all faces by one
                std::transform(_faces.begin(), _faces.end(), _faces.begin(), [](const std::array<size_t, 3> &face) {
                    using namespace util;
                    return face - 1;
                });

            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    bool TetgenAdapter::checkIntegrity(const std::string &filename, char what) const {
        if ((what == 'v' || what == 'a') && _tetgenio.numberofpoints != 0) {
            throw std::runtime_error(
                    "The Polyhedron already has well defined nodes! The information of " + filename
                    + ".node is redundant!");
        } else if ((what == 'f' || what == 'a') && (_tetgenio.numberoftrifaces != 0 || _tetgenio.numberoffacets != 0)) {
            throw std::runtime_error(
                    "The Polyhedron already has well defined faces! The information of " + filename
                    + ".node is redundant!");
        }
        return true;
    }

    void TetgenAdapter::addVertices() {
        _vertices.clear();
        _vertices.reserve(_tetgenio.numberofpoints);
        for (size_t i = 0; i < _tetgenio.numberofpoints * 3; i += 3) {
            _vertices.push_back({_tetgenio.pointlist[i],
                             _tetgenio.pointlist[i + 1],
                             _tetgenio.pointlist[i + 2]});
        }
    }

    void TetgenAdapter::addFacesByTrifaces() {
        _faces.clear();
        _faces.reserve(_tetgenio.numberoftrifaces);
        for (size_t i = 0; i < _tetgenio.numberoftrifaces * 3; i += 3) {
            _faces.push_back({static_cast<size_t>(_tetgenio.trifacelist[i]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 1]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 2])});
        }
    }

    void TetgenAdapter::addFacesByFacetList() {
        _faces.clear();
        _faces.reserve(_tetgenio.numberoffacets);
        for (size_t i = 0; i < _tetgenio.numberoffacets; ++i) {
            _faces.push_back({static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[0]),
                             static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[1]),
                             static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[2])});
        }
    }

}
