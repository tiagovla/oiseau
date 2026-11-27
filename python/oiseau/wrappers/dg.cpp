// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_hexahedron.hpp"
#include "oiseau/dg/nodal/ref_line.hpp"
#include "oiseau/dg/nodal/ref_quadrilateral.hpp"
#include "oiseau/dg/nodal/ref_tetrahedron.hpp"
#include "oiseau/dg/nodal/ref_triangle.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/dg/solver/state.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"

namespace py = nanobind;

namespace oiseau::wrappers {
namespace {
using NumpyCoordinates = py::ndarray<py::numpy, const double, py::c_contig, py::device::cpu>;

xt::xarray<double> copy_coordinates(const NumpyCoordinates& coordinates) {
  const auto* shape = coordinates.shape_ptr();
  auto values = xt::xarray<double>::from_shape(
      xt::xarray<double>::shape_type(shape, shape + coordinates.ndim()));
  std::copy_n(coordinates.data(), coordinates.size(), values.data());
  return values;
}

py::ndarray<py::numpy, double> copy_array(const xt::xarray<double>& values) {
  auto data = std::make_unique_for_overwrite<double[]>(values.size());
  std::copy_n(values.data(), values.size(), data.get());
  py::capsule owner(data.get(),
                    [](void* pointer) noexcept { delete[] static_cast<double*>(pointer); });
  return {data.release(), values.dimension(), values.shape().data(), owner};
}

void bind_nodal_elements(py::module_& m) {
  using namespace oiseau::dg::nodal;

  py::enum_<RefElementType>(m, "RefElementType")
      .value("Line", RefElementType::Line)
      .value("Triangle", RefElementType::Triangle)
      .value("Quadrilateral", RefElementType::Quadrilateral)
      .value("Tetrahedron", RefElementType::Tetrahedron)
      .value("Hexahedron", RefElementType::Hexahedron);

  py::class_<RefElement>(m, "RefElement")
      .def_prop_ro("order", &RefElement::order)
      .def_prop_ro("number_of_nodes", &RefElement::number_of_nodes)
      .def_prop_ro("number_of_face_nodes", &RefElement::number_of_face_nodes)
      .def_prop_ro(
          "r", [](const RefElement& ref) { return copy_array(ref.r()); }, py::rv_policy::reference,
          "Copy of the reference node coordinates.")
      .def_prop_ro(
          "v", [](const RefElement& ref) { return copy_array(ref.v()); }, py::rv_policy::reference,
          "Copy of the Vandermonde matrix at the reference nodes.")
      .def_prop_ro(
          "inv_v", [](const RefElement& ref) { return copy_array(ref.inv_v()); },
          py::rv_policy::reference, "Copy of the inverse Vandermonde matrix.")
      .def_prop_ro(
          "gv", [](const RefElement& ref) { return copy_array(ref.gv()); },
          py::rv_policy::reference,
          "Copy of the gradient Vandermonde matrix at the reference nodes.")
      .def_prop_ro(
          "d", [](const RefElement& ref) { return copy_array(ref.d()); }, py::rv_policy::reference,
          "Copy of the nodal differentiation matrix.")
      .def_prop_ro(
          "mass", [](const RefElement& ref) { return copy_array(ref.mass()); },
          py::rv_policy::reference, "Copy of the reference mass matrix.")
      .def_prop_ro(
          "inv_mass", [](const RefElement& ref) { return copy_array(ref.inv_mass()); },
          py::rv_policy::reference, "Copy of the inverse reference mass matrix.")
      .def(
          "vandermonde",
          [](const RefElement& ref, const NumpyCoordinates& points) {
            const auto& nodes = ref.r();
            if (points.ndim() != nodes.dimension() ||
                (points.ndim() == 2 && points.shape(1) != nodes.shape()[1])) {
              throw py::value_error(
                  "Points must have shape (n,) for line elements or "
                  "(n, reference dimension) for other elements");
            }
            return copy_array(ref.vandermonde(copy_coordinates(points)));
          },
          py::arg("points"), "Evaluate the reference basis at the given points.");

  py::class_<RefLine, RefElement>(m, "RefLine").def(py::init<unsigned>(), py::arg("order"));
  py::class_<RefTriangle, RefElement>(m, "RefTriangle").def(py::init<unsigned>(), py::arg("order"));
  py::class_<RefQuadrilateral, RefElement>(m, "RefQuadrilateral")
      .def(py::init<unsigned>(), py::arg("order"));
  py::class_<RefTetrahedron, RefElement>(m, "RefTetrahedron")
      .def(py::init<unsigned>(), py::arg("order"));
  py::class_<RefHexahedron, RefElement>(m, "RefHexahedron")
      .def(py::init<unsigned>(), py::arg("order"));

  m.def("get_ref_element", &get_ref_element, py::arg("type"), py::arg("order"));
}
}  // namespace

void dg(py::module_& m) {
  using namespace oiseau::dg;
  using solver::LinearDG2D;
  using solver::State;
  using BoundaryCallback =
      std::function<std::vector<double>(const std::vector<double>&, const LinearDG2D::Normal&)>;

  auto nodal = m.def_submodule("nodal", "Nodal Discontinuous Galerkin elements");
  bind_nodal_elements(nodal);

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

  py::class_<Space>(m, "Space")
      .def(py::init<const mesh::Mesh&, unsigned>(), py::arg("mesh"), py::arg("order"),
           py::keep_alive<1, 2>())
      .def(py::init<const mesh::Mesh&, const std::vector<unsigned>&>(), py::arg("mesh"),
           py::arg("order"), py::keep_alive<1, 2>(),
           "Assign one positive polynomial order per top-dimensional cell in mesh order.")
      .def_prop_ro(
          "group_nodes",
          [](const Space& space) {
            py::list nodes;
            for (const auto& group : space.element_groups()) nodes.append(copy_array(group.nodes));
            return nodes;
          },
          "Copies of node coordinates with shape (cells, nodes, dimension), in state group order.");

  py::class_<State>(m, "State")
      .def_prop_ro(
          "groups",
          [](const State& state) {
            py::list groups;
            for (const auto& group : state.groups) groups.append(copy_array(group));
            return groups;
          },
          "Copies of nodal values with shape (cells, nodes, fields), one array per element group.");

  py::class_<LinearDG2D>(m, "LinearDG2D")
      .def(
          "__init__",
          [](LinearDG2D* self, mesh::Mesh& mesh, const Space& space, EquationSystem equations,
             double wave_speed, BoundaryCallback boundary) {
            LinearDG2D::BoundaryState ghost_state;
            if (boundary) {
              ghost_state = [boundary = std::move(boundary)](std::span<const double> interior,
                                                             const LinearDG2D::Normal& normal,
                                                             std::span<double> exterior) {
                const auto values =
                    boundary(std::vector<double>(interior.begin(), interior.end()), normal);
                if (values.size() != exterior.size()) {
                  throw std::invalid_argument(
                      "Boundary state must return one value for every equation field");
                }
                std::ranges::copy(values, exterior.begin());
              };
            }
            new (self) LinearDG2D(
                mesh, space, std::move(equations), [wave_speed](const auto&) { return wave_speed; },
                std::move(ghost_state));
          },
          py::arg("mesh"), py::arg("space"), py::arg("equations"), py::arg("wave_speed"),
          py::arg("boundary_state").none() = py::none(), py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>())
      .def("initial_state", &LinearDG2D::initial_state, py::arg("field"), py::arg("value"));

  m.def(
      "run_rk4",
      [](const LinearDG2D& solver, State state, double initial_time, double time_step,
         std::size_t steps) {
        time::RK4 integrator;
        time::run(integrator, solver, state, initial_time, time_step, steps);
        return state;
      },
      py::arg("solver"), py::arg("state"), py::arg("initial_time"), py::arg("time_step"),
      py::arg("steps"), "Return the state advanced by the given number of fixed RK4 steps.");
}
}  // namespace oiseau::wrappers
