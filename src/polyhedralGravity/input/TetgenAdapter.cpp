#include "TetgenAdapter.h"

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
        _tetgenio.load_node(const_cast<char *>(filename.c_str()));
        _hasNodes = true;
    } else {
        throw std::runtime_error(
                "The Polyhedron already has well defined nodes! The information of " + filename
                + ".node is redundant!");
    }
}

void TetgenAdapter::readFace(const std::string &filename) {
    if (!_hasFaces) {
        _tetgenio.load_face(const_cast<char *>(filename.c_str()));
        _hasFaces = true;
    } else {
        throw std::runtime_error(
                "The Polyhedron already has well defined faces! The information of " + filename
                + ".face is redundant!");
    }
}

void TetgenAdapter::readElements(const std::string &filename) {
    if (!_hasElements) {
        _tetgenio.load_elem(const_cast<char *>(filename.c_str()));
        _hasElements = true;
    } else {
        throw std::runtime_error(
                "The Polyhedron already has well defined elements! The information of " + filename
                + ".ele is redundant!");
    }
}

Polyhedron TetgenAdapter::convertTetgenToPolyhedron() const {
    std::vector<std::array<double, 3>> nodes{};
    nodes.reserve(_tetgenio.numberofpoints);
    for (size_t i = 0; i < _tetgenio.numberofpoints; ++i) {
        nodes.push_back({_tetgenio.pointlist[i],
                         _tetgenio.pointlist[i+1],
                         _tetgenio.pointlist[i+2]});
    }

    std::vector<std::array<size_t, 3>> faces{};
    faces.reserve(_tetgenio.numberoftrifaces);
    for (size_t i = 0; i < _tetgenio.numberoftrifaces; ++i) {
        faces.push_back({static_cast<size_t>(_tetgenio.trifacelist[i]),
                         static_cast<size_t>(_tetgenio.trifacelist[i+1]),
                         static_cast<size_t>(_tetgenio.trifacelist[i+2])});
    }
    return {nodes, faces};
}


