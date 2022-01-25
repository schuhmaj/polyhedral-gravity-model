#include <iostream>

#include "spdlog/spdlog.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/Gravity.h"

int main(int argc, char *argv[]) {
    SPDLOG_INFO("The answer to your question is 42!");

    std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("../example-config/example.yaml");
    auto poly = config->getDataSource()->getPolyhedron();

    Gravity grav{poly};

    grav.calculate();

    SPDLOG_INFO("Finished.");
    return 0;
}