// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <pybind11/cast.h>

#include <iostream>
#include <ostream>
#include <ranges>

#include "fmt/base.h"
#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/plotting/triplot.hpp"

using namespace oiseau::io;
using namespace oiseau::mesh;

namespace py = pybind11;
using namespace pybind11::literals;

void two_dimensional() {
  Mesh mesh = gmsh_read_from_path("demo/meshes/mesh2d_split_triangles_quads.msh");
  plt::scoped_interpreter guard{};

  auto [fig, ax] = plt::subplots(1, 1, py::arg("constrained_layout") = true,
                                 py::arg("figsize") = py::make_tuple(8.0, 8.0));

  oiseau::plotting::triplot(ax, mesh);

  plt::show();
}

void three_dimensional() {
  Mesh mesh = gmsh_read_from_path("demo/meshes/mesh3d_tetra_hexa_block.msh");
  auto& topology = mesh.topology();

  auto cell_types = topology.cell_types();
  auto conn = topology.conn();

  for (auto [ct, vertices] : std::views::zip(cell_types, conn)) {
    std::cout << ct->name() << std::endl;
    fmt::print("{}", vertices);
  }
}

int main() {
  two_dimensional();
  three_dimensional();
}
