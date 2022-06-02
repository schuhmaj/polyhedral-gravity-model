#include <chrono>
#include "polyhedralGravity/input/ConfigSource.h"
#include "polyhedralGravity/input/YAMLConfigReader.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/output/Logging.h"

int main(int argc, char *argv[]) {
    using namespace polyhedralGravity;
    if (argc != 2) {
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "Wrong program call! "
                                                                                 "Please use the program like that:\n"
                                                                                 "./polyhedralGravity [YAML-Configuration-File]\n");
        return 0;
    }

    try {

        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>(argv[1]);
        auto poly = config->getDataSource()->getPolyhedron();
        auto density = config->getDensity();
        auto point = config->getPointsOfInterest()[0];

        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "The calculation started...");
        auto start = std::chrono::high_resolution_clock::now();
        auto result = GravityModel::evaluate(poly, density, point);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "The calculation finished. It took {} microseconds.", ms.count());

        //The results
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "V= {}", result.gravitationalPotential);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vx= {}", result.acceleration[0]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vy= {}", result.acceleration[1]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vz= {}", result.acceleration[2]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vxx= {}", result.gradiometricTensor[0]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vyy= {}", result.gradiometricTensor[1]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vzz= {}", result.gradiometricTensor[2]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vxy= {}", result.gradiometricTensor[3]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vxz= {}", result.gradiometricTensor[4]);
        SPDLOG_LOGGER_INFO(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(),
                           "Vyz= {}", result.gradiometricTensor[5]);
        return 0;

    } catch (const std::exception &e) {
        SPDLOG_LOGGER_ERROR(PolyhedralGravityLogger::DEFAULT_LOGGER.getLogger(), "{}", e.what());
        return -1;
    }
}