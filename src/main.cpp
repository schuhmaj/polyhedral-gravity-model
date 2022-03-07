#include <chrono>
#include "spdlog/spdlog.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/Gravity.h"

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;
    SPDLOG_INFO("The answer to your question is 42!");

    std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("../example-config/example.yaml");
    auto poly = config->getDataSource()->getPolyhedron();
    auto density = config->getDensity();

    Gravity grav{poly, density};

    auto start = std::chrono::high_resolution_clock::now();
    grav.calculate();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    SPDLOG_INFO("The calculation took {} microseconds", ms.count());

    SPDLOG_INFO("Finished.");
    return 0;
}