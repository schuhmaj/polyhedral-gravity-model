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
        return convertTetgenToPolyhedron();
    }

    void TetgenAdapter::readNode(const std::string &filename) {
        if (this->checkIntegrity(filename, 'v')) {
            try {
                _tetgenio.load_node(const_cast<char *>(filename.c_str()));
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
            } catch (...) {
                throw std::runtime_error(
                        "The faces were not read because of an error in Tetgen! This could indicate several "
                        "issues, e. g. issues with the node assignment like they appear if either no nodes were "
                        "read in at all or if no assignment was possible."
                );
            }
        }
    }

    void TetgenAdapter::readElements(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_elem(const_cast<char *>(filename.c_str()));
            } catch (...) {
                throw std::runtime_error(
                        "The elements were not read because of an error in Tetgen!"
                );
            }
        }
    }

    void TetgenAdapter::readOff(const std::string &filename) {
        if (this->checkIntegrity(filename, 'a')) {
            try {
                _tetgenio.load_off(const_cast<char *>(filename.c_str()));
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

    Polyhedron TetgenAdapter::convertTetgenToPolyhedron() const {
        std::vector<std::array<double, 3>> nodes{};
        nodes.reserve(_tetgenio.numberofpoints);
        for (size_t i = 0; i < _tetgenio.numberofpoints * 3; i += 3) {
            nodes.push_back({_tetgenio.pointlist[i],
                             _tetgenio.pointlist[i + 1],
                             _tetgenio.pointlist[i + 2]});
        }

        //Trifacet versio
        std::vector<std::array<size_t, 3>> faces{};
        faces.reserve(_tetgenio.numberoftrifaces);
        for (size_t i = 0; i < _tetgenio.numberoftrifaces * 3; i += 3) {
            faces.push_back({static_cast<size_t>(_tetgenio.trifacelist[i]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 1]),
                             static_cast<size_t>(_tetgenio.trifacelist[i + 2])});
        }
        //Polygon version
        faces.reserve(_tetgenio.numberoffacets);
        for (size_t i = 0; i < _tetgenio.numberoffacets; i += 3) {
            //TODO Error Handling if polygonlist > 1
            faces.push_back({static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[0]),
                             static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[1]),
                             static_cast<size_t>(_tetgenio.facetlist[i].polygonlist->vertexlist[2])});
        }
        return {nodes, faces};
    }

}
