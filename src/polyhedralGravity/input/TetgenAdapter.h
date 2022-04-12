#pragma once

#include "DataSource.h"
#include "tetgen.h"
#include <map>
#include <string>
#include <functional>
#include <utility>
#include <exception>
#include <stdexcept>

namespace polyhedralGravity {

/**
 * An adapter to the Tetgen Library. This is the interface between tetgen's data structures and I/O
 * capabilities and the here implemented polyhedral gravity model.
 * The adapter always converts tetgen's datastructure into the structure utilized by the Polyhedral Gravity Model.
 *
 * The Adapter further keeps en eye on the already read in files in order to give feedback if data is in conflict.
 */
    class TetgenAdapter : public DataSource {

        /**
         * Delegates the call to the library tetgen
         */
        tetgenio _tetgenio;

        /**
         * The file names from which to read in the polyhedron
         */
        const std::vector<std::string> _fileNames;

    public:

        explicit TetgenAdapter(std::vector<std::string> fileNames)
                : _tetgenio(),
                  _fileNames(std::move(fileNames)) {}

        /**
         * Use this function to get a Polyhedron.
         * This functions consists of two steps. First, the Adapter will delegate I/O to the tetgen library and
         * read in the Polyhedron data in the library's datastructure. Second, tetgen's datastructure is then
         * converted to a Polyhedron.
         * @return a Polyhedron
         */
        Polyhedron getPolyhedron() override;

        /**
         * Reads nodes from a .node file
         * @param filename of the input source without suffix
         * @throws an exception if the nodes already have been defined
         */
        void readNode(const std::string &filename);

        /**
         * Reads faces from a .face file
         * @param filename of the input source without suffix
         * @throws an exception if the faces already have been defined
         */
        void readFace(const std::string &filename);

        /**
         * Reads elements from a .ele file
         * @param filename of the input source without suffix
         * @throws an exception if the elements already have been defined
         */
        void readElements(const std::string &filename);


        /**
         * Reads elements from a .off file (Geomview’s polyhedral file format)
         * @param filename of the input source without suffix
         * @throws an exception if the elements already have been defined
         */
        void readOff(const std::string &filename);


        /**
         * Reads elements from a .ply file (Polyhedral file format)
         * @param filename of the input source without suffix
         * @throws an exception if the elements already have been defined
         */
        void readPly(const std::string &filename);


        /**
         * Reads elements from a .stl file (Stereolithography format)
         * @param filename of the input source without suffix
         * @throws an exception if the elements already have been defined
         */
        void readStl(const std::string &filename);

        /**
         * Reads elements from a .mesh file (Medit’s mesh file format)
         * @param filename of the input source without suffix
         * @throws an exception if the elements already have been defined
         */
        void readMesh(const std::string &filename);



    private:

        /**
         * The map contains the file suffix and the corresponding operation applicable for this suffix.
         * File suffix --> Operation on this instance
         */
        const std::map<std::string, std::function<void(const std::string &name)>> _suffixToOperation{
                {"node", [this](const std::string &name) { this->readNode(name); }},
                {"face", [this](const std::string &name) { this->readFace(name); }},
                {"ele",  [this](const std::string &name) { this->readElements(name); }},
                {"off",  [this](const std::string &name) { this->readOff(name); }},
                {"ply",  [this](const std::string &name) { this->readPly(name); }},
                {"stl",  [this](const std::string &name) { this->readStl(name); }},
                {"mesh",  [this](const std::string &name) { this->readMesh(name); }}
        };

        /**
         * Checks if the polyhedron is integer and not already defined by other properties
         * @param filename - string with the current read file, for more detailed exceptions
         * @return true if everything is ok
         * @throws an exception if not
         */
        bool checkIntegrity(const std::string &filename) const;

        /**
         * Converts the tetgenio structure to a more slim Polyhedron used by this implementation.
         * @return a Polyhedron
         */
        [[nodiscard]] Polyhedron convertTetgenToPolyhedron() const;
    };

}
