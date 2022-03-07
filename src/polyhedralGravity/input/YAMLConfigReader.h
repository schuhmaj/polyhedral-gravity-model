#pragma once

#include "ConfigSource.h"
#include "TetgenAdapter.h"
#include "yaml-cpp/yaml.h"
#include <string>
#include <vector>
#include <array>

namespace polyhedralGravity {

    class YAMLConfigReader : public ConfigSource {

        /*
         * The following static variables contain the names of the YAML nodes.
         * Note: C++20 would allow constexpr std::string which would be more appropriate instead of char[]
         */
        static constexpr char ROOT[] = "gravityModel";
        static constexpr char INPUT[] = "input";
        static constexpr char INPUT_POLYHEDRON[] = "polyhedron";
        static constexpr char INPUT_DENSITY[] = "density";
        static constexpr char INPUT_POINTS[] = "points";

        const YAML::Node _file;

    public:

        /**
         * Creates a new YAML Config Reader.
         * @param filename - a reference to a string
         * @throws an exception if the file is malformed or cannot be loaded or if the ROOT node is not found
         */
        explicit YAMLConfigReader(const std::string &filename)
                : _file{YAML::LoadFile(filename)} {
            if (!_file[ROOT]) {
                throw std::runtime_error{"The YAML file does not contain a specification for the \"gravityModel\"!"};
            }
        }

        /**
         * Creates a new YAML Config Reader.
         * @param filename - a movable string
         * @throws an exception if the file is malformed or cannot be loaded or if the ROOT node is not found
         */
        explicit YAMLConfigReader(std::string &&filename)
                : _file{YAML::LoadFile(filename)} {
            if (!_file[ROOT]) {
                throw std::runtime_error{"The YAML file does not contain a specification for the \"gravityModel\"!"};
            }
        }

        double getDensity() override;

        std::vector<std::array<double, 3>> getPointsOfInterest() override;

        std::shared_ptr<DataSource> getDataSource() override;


    };

}
