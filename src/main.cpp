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

    std::array<size_t, 1000> time{};
    for (int i = 0; i < 1000; ++i) {
        //SPDLOG_INFO("The calculation started.");
        auto start = std::chrono::high_resolution_clock::now();
        auto result = GravityModel::evaluate(poly, density, point);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end - start;
        auto ms = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        time[i] = ms.count();
    }

    double avg = 0;
    for (int i = 0; i < 1000; ++i) {
        avg += time[i];
    }
    avg /= 1000.0;
    SPDLOG_INFO("Average {} microseconds.", avg);

    return 0;
}
