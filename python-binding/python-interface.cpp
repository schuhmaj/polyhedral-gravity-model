#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityModel.h"


namespace py = pybind11;

//TODO Wrapper required! Otherwise segfault if parallelized!
polyhedralGravity::GravityModelResult eval(const polyhedralGravity::Polyhedron &polyhedron) {
    return polyhedralGravity::GravityModel::evaluate(polyhedron, 2760.0);
}

PYBIND11_MODULE(polyhedral_gravity, m) {
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron at a given computation point P";

    py::class_<polyhedralGravity::GravityModelResult>(m, "GravityModelResult")
            .def(py::init<>())
            .def(py::init<double, const std::array<double, 3> &, const std::array<double, 6> &>())
            .def_readwrite("potential", &polyhedralGravity::GravityModelResult::gravitationalPotential)
            .def_readwrite("acceleration", &polyhedralGravity::GravityModelResult::gravitationalPotentialDerivative)
            .def_readwrite("tensor", &polyhedralGravity::GravityModelResult::gradiometricTensor);

    py::class_<polyhedralGravity::Polyhedron>(m, "Polyhedron")
            .def(py::init<>())
            .def(py::init<std::vector<std::array<double, 3>>, std::vector<std::array<size_t, 3>>>())
            .def("vertices", &polyhedralGravity::Polyhedron::getVertices)
            .def("faces", &polyhedralGravity::Polyhedron::getFaces);

    m.def("evaluate", &eval,
          "Evaluate the full gravity tensor for a given constant density polyhedron at a given computation point P",
          py::arg("polyhedron"));

}