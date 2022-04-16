#include <chrono>
#include "spdlog/spdlog.h"
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/GravityModel.h"

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;
    if (argc != 2) {
        SPDLOG_INFO("Wrong program call! Please use the program like that:\n"
                    "./polyhedralGravity [YAML-Configuration-File]\n");
        return 0;
    }

    std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
    auto poly = config->getDataSource()->getPolyhedron();
    auto density = config->getDensity();

    SPDLOG_INFO("The calculation started.");
    auto start = std::chrono::high_resolution_clock::now();
    GravityModel::evaluate(poly, density, {0.0, 0.0, 0.0});
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    SPDLOG_INFO("The calculation finished. It took {} microseconds.", ms.count());
    return 0;
}