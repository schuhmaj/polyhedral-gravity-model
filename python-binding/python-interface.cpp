#include <tuple>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityModel.h"
#include "polyhedralGravity/input/TetgenAdapter.h"


namespace py = pybind11;

std::tuple<double, std::array<double, 3>, std::array<double, 6>> convertToTuple(
        const polyhedralGravity::GravityModelResult &result) {
    return std::make_tuple(result.gravitationalPotential,
                           result.gravitationalPotentialDerivative,
                           result.gradiometricTensor);
}

std::vector<std::tuple<double, std::array<double, 3>, std::array<double, 6>>> convertToTupleVector(
        const std::vector<polyhedralGravity::GravityModelResult> &resultVector) {
    std::vector<std::tuple<double, std::array<double, 3>, std::array<double, 6>>> resVectorTuple{resultVector.size()};
    std::transform(resultVector.cbegin(), resultVector.cend(), resVectorTuple.begin(), &convertToTuple);
    return resVectorTuple;
}

PYBIND11_MODULE(polyhedral_gravity, m) {
    using namespace polyhedralGravity;
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron which consists of some "
              "vertices and triangular faces at a given computation point P";

    /*
     * Methods for vertices and faces from the script context
     */

    m.def("evaluate",
          [](const std::vector<std::array<double, 3>> &vertices, const std::vector<std::array<size_t, 3>> &faces,
             double density, const std::array<double, 3> &computationPoint) {
              return convertToTuple(GravityModel::evaluate({vertices, faces}, density, computationPoint));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at a given computation point P",
          py::arg("vertices"), py::arg("faces"), py::arg("density"), py::arg("computation_point"));

    m.def("evaluate",
          [](const std::vector<std::array<double, 3>> &vertices, const std::vector<std::array<size_t, 3>> &faces,
             double density, const std::vector<std::array<double, 3>> &computationPoints) {
              return convertToTupleVector(GravityModel::evaluate({vertices, faces}, density, computationPoints));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at multiple given computation points",
          py::arg("vertices"), py::arg("faces"), py::arg("density"), py::arg("computation_points"));

    /*
     * Methods for vertices and faces from .node and .face files via TetGen
     */

    m.def("evaluate",
          [](const std::string &node_file, const std::string &face_file,
             double density, const std::array<double, 3> &computationPoint) {
              TetgenAdapter tetgen{{node_file, face_file}};
              return convertToTuple(GravityModel::evaluate(tetgen.getPolyhedron(), density, computationPoint));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices"
          "and triangular faces at a given computation point P. The vertices and faces are read from .node and .face "
          "files.",
          py::arg("node_file"), py::arg("face_file"), py::arg("density"), py::arg("computation_point"));

    m.def("evaluate",
          [](const std::string &node_file, const std::string &face_file,
             double density, const std::vector<std::array<double, 3>> &computationPoints) {
              TetgenAdapter tetgen{{node_file, face_file}};
              return convertToTupleVector(GravityModel::evaluate(tetgen.getPolyhedron(), density, computationPoints));
          },
          "Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices "
          "and triangular faces at multiple given computation points. The vertices and faces are read from "
          ".node and .face files.",
          py::arg("node_file"), py::arg("face_file"), py::arg("density"), py::arg("computation_points"));
}
