#include "TetgenAdapter.h"

Polyhedron TetgenAdapter::getPolyhedron() {
    //1. Step: Read in from files
    for (const auto &fileName : _fileNames) {
        size_t pos = fileName.find_last_of('.');
        std::string name = fileName.substr(0, pos);
        std::string suffix = fileName.substr(pos + 1);
        _suffixToOperation.at(suffix)(name);
    }

    //2. Convert to tetgenio to Polyhedron

    //TODO Conversion to Polyhedron
    return Polyhedron{};
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


