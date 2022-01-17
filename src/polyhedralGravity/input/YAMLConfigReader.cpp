#include "YAMLConfigReader.h"

std::vector<std::array<double, 3>> YAMLConfigReader::getPointsOfInterest() {
    if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POINTS]) {
        return _file[ROOT][INPUT][INPUT_POINTS].as<std::vector<std::array<double, 3>>>();
    } else {
        throw std::runtime_error{"There happened an error parsing the points of interest from the YAML config file!"};
    }
}

std::shared_ptr<DataSource> YAMLConfigReader::getDataSource() {
    if (_file[ROOT][INPUT] && _file[ROOT][INPUT][INPUT_POLYHEDRON]) {
        auto vectorOfFiles = _file[ROOT][INPUT][INPUT_POLYHEDRON].as<std::vector<std::string>>();
        return std::make_shared<TetgenAdapter>(vectorOfFiles);
    } else {
        throw std::runtime_error{"There happened an error parsing the DataSource of the Polyhedron from the config file"};
    }
}
