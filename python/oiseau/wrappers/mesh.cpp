// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <vector>

#include "oiseau/mesh/cell.hpp"

namespace py = pybind11;

namespace oiseau::wrappers {
void mesh(py::module &m) {
  py::class_<mesh::Cell, std::unique_ptr<mesh::Cell, py::nodelete>> cell(m, "Cell");
  cell.def("facet", &mesh::Cell::facet);
  cell.def("name", &mesh::Cell::name);
  cell.def(
      "get_entity_vertices",
      [](mesh::Cell &self, int dim) {
        std::vector<std::vector<int>> cpp_result = self.get_entity_vertices(dim);
        return py::array(py::cast(cpp_result));
      },
      py::arg("dim"), "Get entity vertices.");
  m.def("get_cell_type", &mesh::get_cell_type);
}

}  // namespace oiseau::wrappers
