#include <iostream>

#include "spdlog/spdlog.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"

int main(int argc, char *argv[]) {
    SPDLOG_INFO("The answer to your question is 42!");

    std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("../example-config/example.yaml");
    auto points = config->getDataSource()->getPolyhedron();

    return 0;
}