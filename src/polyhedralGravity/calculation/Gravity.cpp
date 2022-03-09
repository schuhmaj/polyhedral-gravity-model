#include "Gravity.h"

namespace polyhedralGravity {

    void Gravity::calculate() {
        SPDLOG_INFO("Calculate...");
        auto x = calculateGij();
    }

    std::vector<std::array<std::array<double, 3>, 3>> Gravity::calculateGij() {
        std::vector<std::array<std::array<double, 3>, 3>> g{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), g.begin(),
                       [&](const auto &face) -> std::array<std::array<double, 3>, 3> {
                           using util::operator-;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);
                           return {node1 - node0, node2 - node1, node0 - node2};
                       });
        return g;
    }

    std::vector<std::array<double, 3>>
    Gravity::calculatePlaneUnitNormals(const std::vector<std::array<std::array<double, 3>, 3>> &g) {
        std::vector<std::array<double, 3>> planeUnitNormal{g.size()};
        std::transform(g.cbegin(), g.cend(), planeUnitNormal.begin(), [](const auto &gi) -> std::array<double, 3> {
            using namespace util;
            const std::array<double, 3> crossProduct = cross(gi[0], gi[1]);
            const double norm = euclideanNorm(crossProduct);
            return crossProduct / norm;
        });
        return planeUnitNormal;
    }

    std::vector<std::array<std::array<double, 3>, 3>>
    Gravity::calculateSegmentUnitNormals(const std::vector<std::array<std::array<double, 3>, 3>> &g,
                                         const std::vector<std::array<double, 3>> &planeUnitNormals) {
        std::vector<std::array<std::array<double, 3>, 3>> segmentUnitNormal{g.size()};
        //Outer "loop" over G_i (running i) and N_i calculating n_i
        std::transform(g.cbegin(), g.cend(), planeUnitNormals.cbegin(), segmentUnitNormal.begin(),
                       [](const std::array<std::array<double, 3>, 3> &gi, const std::array<double, 3> &Ni)
                               -> std::array<std::array<double, 3>, 3> {
                           std::array<std::array<double, 3>, 3> ni{};
                           //Inner "loop" over G_ij (fixed i, running j) with parameter N_i calculating n_ij
                           std::transform(gi.cbegin(), gi.end(), ni.begin(),
                                          [&Ni](const auto &gij) -> std::array<double, 3> {
                                              using namespace util;
                                              const std::array<double, 3> crossProduct = cross(gij, Ni);
                                              const double norm = euclideanNorm(crossProduct);
                                              return crossProduct / norm;
                                          });
                           return ni;
                       });
        return segmentUnitNormal;
    }

    std::vector<double> Gravity::calculateSigmaP(const std::vector<std::array<double, 3>> &planeUnitNormals) {
        std::vector<double> sigmaP(planeUnitNormals.size(), 0.0);
        std::transform(planeUnitNormals.cbegin(), planeUnitNormals.cend(), _polyhedron.getFaces().begin(), sigmaP.begin(),
                       [&](const std::array<double, 3> &ni, const std::array<size_t, 3> &gi) {
                           using namespace util;
                           const auto &node = _polyhedron.getNode(gi[0]);
                           return sgn(dot(ni, node * -1.0)) * -1.0;
                       });
        return sigmaP;
    }


    std::vector<HessianPlane> Gravity::calculateFaceToHessianPlane(const std::array<double, 3> &p) {
        std::vector<HessianPlane> hessianPlane{_polyhedron.countFaces()};
        std::transform(_polyhedron.getFaces().cbegin(), _polyhedron.getFaces().cend(), hessianPlane.begin(),
                       [&](const auto &face) -> HessianPlane {
                           using namespace util;
                           const auto &node0 = _polyhedron.getNode(face[0]);
                           const auto &node1 = _polyhedron.getNode(face[1]);
                           const auto &node2 = _polyhedron.getNode(face[2]);

                           return computeHessianPlane(node0, node1, node2, p);
                       });
        return hessianPlane;
    }

    HessianPlane Gravity::computeHessianPlane(const std::array<double, 3> &p, const std::array<double, 3> &q,
                                              const std::array<double, 3> &r, const std::array<double, 3> &origin) {
        using namespace util;
        const auto crossProduct = cross(p - q, p - r);
        const auto res = (origin - p) * crossProduct;
        const double d = res[0] + res[1] + res[2];

        return {crossProduct[0], crossProduct[1], crossProduct[2], d};
    }

    std::vector<double> Gravity::calculatePlaneDistance(const std::vector<HessianPlane> &plane) {
        std::vector<double> planeDistances(plane.size(), 0.0);
        std::transform(plane.cbegin(), plane.cend(), planeDistances.begin(),
                       [](const HessianPlane &plane) -> double {
                           return std::abs(
                                   plane.d / std::sqrt(plane.a * plane.a + plane.b * plane.b + plane.c * plane.c));
                       });
        return planeDistances;
    }


}
