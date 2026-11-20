// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/mesh.hpp"

#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "oiseau/mesh/cell.hpp"

namespace py = pybind11;

namespace oiseau::wrappers {
void mesh(py::module& m) {
  using mesh::CellKind;
  py::enum_<CellKind>(m, "CellKind")
      .value("Point", CellKind::Point)
      .value("Interval", CellKind::Interval)
      .value("Triangle", CellKind::Triangle)
      .value("Quadrilateral", CellKind::Quadrilateral)
      .value("Tetrahedron", CellKind::Tetrahedron)
      .value("Hexahedron", CellKind::Hexahedron);

  py::class_<mesh::Cell, std::unique_ptr<mesh::Cell, py::nodelete>> cell(m, "Cell");
  cell.def("facet", &mesh::Cell::facet, py::return_value_policy::reference);
  cell.def("name", &mesh::Cell::name);
  cell.def(
      "get_entity_vertices",
      [](mesh::Cell& self, int dim) {
        if (dim < 0 || dim > self.dimension()) throw py::value_error("Invalid entity dimension");
        std::vector<std::vector<int>> cpp_result = self.get_entity_vertices(dim);
        return py::array(py::cast(cpp_result));
      },
      py::arg("dim"), "Get entity vertices.");
  m.def("get_cell_type", &mesh::get_cell_type, py::return_value_policy::reference);

  py::class_<mesh::Mesh>(m, "Mesh").def(
      py::init([](const py::array_t<double, py::array::c_style | py::array::forcecast>& points,
                  std::vector<std::vector<std::size_t>> cells, const std::vector<CellKind>& kinds) {
        if (points.ndim() != 2 || points.shape(1) < 1 || points.shape(1) > 3) {
          throw py::value_error(
              "Points must have shape (vertex count, dimension) with 1 to 3 coordinates");
        }
        std::vector<mesh::CellType> types;
        types.reserve(kinds.size());
        for (const auto kind : kinds) types.push_back(mesh::get_cell_type(kind));
        std::vector<double> coordinates(points.data(), points.data() + points.size());
        return mesh::Mesh(mesh::Topology(std::move(cells), std::move(types)),
                          mesh::Geometry(std::move(coordinates), points.shape(1)));
      }),
      py::arg("points"), py::arg("cells"), py::arg("cell_kinds"));
}

}  // namespace oiseau::wrappers
