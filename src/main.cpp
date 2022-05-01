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
    auto point = config->getPointsOfInterest()[0];

    SPDLOG_INFO("The calculation started.");
    auto start = std::chrono::high_resolution_clock::now();
    auto result = GravityModel::evaluate(poly, density, point);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    SPDLOG_INFO("The calculation finished. It took {} microseconds.", ms.count());

    //The results
    SPDLOG_INFO("V= {}", result.gravitationalPotential);
    SPDLOG_INFO("Vx= {}", result.gravitationalPotentialDerivative[0]);
    SPDLOG_INFO("Vy= {}", result.gravitationalPotentialDerivative[1]);
    SPDLOG_INFO("Vz= {}", result.gravitationalPotentialDerivative[2]);
    SPDLOG_INFO("Vxx= {}", result.gradiometricTensor[0]);
    SPDLOG_INFO("Vyy= {}", result.gradiometricTensor[1]);
    SPDLOG_INFO("Vzz= {}", result.gradiometricTensor[2]);
    SPDLOG_INFO("Vxy= {}", result.gradiometricTensor[3]);
    SPDLOG_INFO("Vxz= {}", result.gradiometricTensor[4]);
    SPDLOG_INFO("Vyz= {}", result.gradiometricTensor[5]);
    return 0;
}