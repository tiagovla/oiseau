// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"

namespace py = pybind11;

namespace oiseau::wrappers {
namespace {
py::array_t<double> copy_array(const xt::xarray<double>& values) {
  py::array_t<double> result(values.shape());
  std::copy_n(values.data(), values.size(), result.mutable_data());
  return result;
}
}  // namespace

void dg(py::module& m) {
  using namespace oiseau::dg;
  using solver::DGState;
  using solver::LinearDG2D;
  using BoundaryCallback =
      std::function<std::vector<double>(const std::vector<double>&, const LinearDG2D::Normal&)>;

  py::class_<Field>(m, "Field").def(py::init<std::string>(), py::arg("name"));
  py::class_<Equation>(m, "Equation");
  py::class_<LinearExpression>(m, "LinearExpression")
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(-py::self)
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self == py::self);
  m.def("dt", [](const Field& field) { return LinearExpression(dt(field)); });
  m.def("dx", [](const Field& field) { return LinearExpression(dx(field)); });
  m.def("dy", [](const Field& field) { return LinearExpression(dy(field)); });

  py::class_<EquationSystem>(m, "EquationSystem")
      .def(py::init<std::vector<Field>, std::vector<Equation>>(), py::arg("fields"),
           py::arg("equations"))
      .def("field_index", &EquationSystem::field_index, py::arg("field"));

  py::class_<DGSpace>(m, "DGSpace")
      .def(py::init<const mesh::Mesh&, unsigned>(), py::arg("mesh"), py::arg("order"),
           py::keep_alive<1, 2>())
      .def(py::init<const mesh::Mesh&, const std::vector<unsigned>&>(), py::arg("mesh"),
           py::arg("order"), py::keep_alive<1, 2>(),
           "Assign one positive polynomial order per top-dimensional cell in mesh order.")
      .def_property_readonly(
          "group_nodes",
          [](const DGSpace& space) {
            py::list nodes;
            for (const auto& group : space.element_groups()) nodes.append(copy_array(group.nodes));
            return nodes;
          },
          "Copies of node coordinates with shape (cells, nodes, dimension), in state group order.");

  py::class_<DGState>(m, "DGState")
      .def_property_readonly(
          "groups",
          [](const DGState& state) {
            py::list groups;
            for (const auto& group : state.groups) groups.append(copy_array(group));
            return groups;
          },
          "Copies of nodal values with shape (cells, nodes, fields), one array per element group.");

  py::class_<LinearDG2D>(m, "LinearDG2D")
      .def(py::init([](mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
                       double wave_speed, BoundaryCallback boundary) {
             LinearDG2D::BoundaryState ghost_state;
             if (boundary) {
               ghost_state = [boundary = std::move(boundary)](std::span<const double> interior,
                                                              const LinearDG2D::Normal& normal) {
                 return boundary(std::vector<double>(interior.begin(), interior.end()), normal);
               };
             }
             return LinearDG2D(
                 mesh, space, std::move(equations),
                 [wave_speed](const auto&) { return wave_speed; }, std::move(ghost_state));
           }),
           py::arg("mesh"), py::arg("space"), py::arg("equations"), py::arg("wave_speed"),
           py::arg("boundary_state") = py::none(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
      .def("initial_state", &LinearDG2D::initial_state, py::arg("field"), py::arg("value"));

  m.def(
      "run_rk4",
      [](const LinearDG2D& solver, DGState state, double initial_time, double time_step,
         std::size_t steps) {
        time::RK4 integrator;
        time::run(integrator, solver, state, initial_time, time_step, steps);
        return state;
      },
      py::arg("solver"), py::arg("state"), py::arg("initial_time"), py::arg("time_step"),
      py::arg("steps"), "Return the state advanced by the given number of fixed RK4 steps.");
}
}  // namespace oiseau::wrappers
