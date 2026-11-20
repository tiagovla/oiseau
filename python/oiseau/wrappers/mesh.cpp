// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/mesh.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <new>
#include <utility>
#include <vector>

#include "oiseau/mesh/cell.hpp"

namespace py = nanobind;

namespace oiseau::wrappers {
void mesh(py::module_& m) {
  using mesh::CellKind;
  auto cell_kind_value = [](CellKind kind) { return static_cast<int>(kind); };
  py::enum_<CellKind>(m, "CellKind")
      .value("Point", CellKind::Point)
      .value("Interval", CellKind::Interval)
      .value("Triangle", CellKind::Triangle)
      .value("Quadrilateral", CellKind::Quadrilateral)
      .value("Tetrahedron", CellKind::Tetrahedron)
      .value("Hexahedron", CellKind::Hexahedron)
      .def("__int__", cell_kind_value)
      .def("__index__", cell_kind_value)
      .def("__hash__", cell_kind_value);

  py::class_<mesh::Cell> cell(m, "Cell", py::never_destruct());
  cell.def("facet", &mesh::Cell::facet, py::rv_policy::reference);
  cell.def("name", &mesh::Cell::name);
  cell.def(
      "get_entity_vertices",
      [](mesh::Cell& self, int dim) {
        if (dim < 0 || dim > self.dimension()) throw py::value_error("Invalid entity dimension");
        return py::module_::import_("numpy").attr("asarray")(
            py::cast(self.get_entity_vertices(dim)));
      },
      py::arg("dim"), "Get entity vertices.");
  m.def("get_cell_type", &mesh::get_cell_type, py::rv_policy::reference);

  py::class_<mesh::Mesh>(m, "Mesh").def(
      "__init__",
      [](mesh::Mesh* self, py::handle points, std::vector<std::vector<std::size_t>> cells,
         const std::vector<CellKind>& kinds) {
        auto array = py::module_::import_("numpy").attr("asarray")(
            points, py::arg("dtype") = "float64", py::arg("order") = "C");
        auto vertices = py::cast<py::ndarray<const double, py::c_contig>>(array);
        if (vertices.ndim() != 2 || vertices.shape(1) < 1 || vertices.shape(1) > 3) {
          throw py::value_error(
              "Points must have shape (vertex count, dimension) with 1 to 3 coordinates");
        }
        std::vector<mesh::CellType> types;
        types.reserve(kinds.size());
        for (const auto kind : kinds) types.push_back(mesh::get_cell_type(kind));
        std::vector<double> coordinates(vertices.data(), vertices.data() + vertices.size());
        new (self) mesh::Mesh(mesh::Topology(std::move(cells), std::move(types)),
                              mesh::Geometry(std::move(coordinates), vertices.shape(1)));
      },
      py::arg("points"), py::arg("cells"), py::arg("cell_kinds"));
}

}  // namespace oiseau::wrappers
