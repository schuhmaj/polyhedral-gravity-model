#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "polyhedralGravity/model/Polyhedron.h"
#include "polyhedralGravity/model/GravityModelData.h"
#include "polyhedralGravity/calculation/GravityModel.h"


namespace py = pybind11;

PYBIND11_MODULE(polyhedral_gravity, m) {
    using namespace polyhedralGravity;
    m.doc() = "Computes the full gravity tensor for a given constant density polyhedron at a given computation point P";

    py::class_<GravityModelResult>(m, "GravityModelResult")
            .def(py::init<>())
            .def(py::init<double, const std::array<double, 3> &, const std::array<double, 6> &>())
            .def_readwrite("potential", &GravityModelResult::gravitationalPotential)
            .def_readwrite("acceleration", &GravityModelResult::gravitationalPotentialDerivative)
            .def_readwrite("tensor", &GravityModelResult::gradiometricTensor);

    py::class_<polyhedralGravity::Polyhedron>(m, "Polyhedron")
            .def(py::init<>())
            .def(py::init<std::vector<std::array<double, 3>>, std::vector<std::array<size_t, 3>>>())
            .def("vertices", &Polyhedron::getVertices)
            .def("faces", &Polyhedron::getFaces);

    m.def("evaluate", py::overload_cast<const Polyhedron &, double, const Array3 &>(&GravityModel::evaluate),
          "Evaluate the full gravity tensor for a given constant density polyhedron "
          "at a given computation point P",
          py::arg("polyhedron"), py::arg("density"), py::arg("computation_point"));

    m.def("evaluate",
          py::overload_cast<const Polyhedron &, double, const std::vector<Array3> &>(&GravityModel::evaluate),
          "Evaluate the full gravity tensor for a given constant density polyhedron "
          "at multiple given computation points",
          py::arg("polyhedron"), py::arg("density"), py::arg("vector_computation_points"));
}
