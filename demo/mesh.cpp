// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/mesh.hpp"

#include <fmt/base.h>
#include <fmt/ranges.h>
#include <pybind11/cast.h>

#include <cstddef>
#include <utility>
#include <vector>
#include <xtensor/containers/xbuffer_adaptor.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/plotting/triplot.hpp"
#include "xtensor/containers/xadapt.hpp"

using namespace oiseau::mesh;
int main() {
  std::vector<double> x = {-2, -1, 0, -1, 0, 1, -2, 1, 2, -1, 2, 1};
  std::vector<std::vector<std::size_t>> conn{{0, 1, 2}, {0, 2, 3}, {1, 4, 5, 2}};

  auto triangle_type = get_cell_type(CellKind::Triangle);
  auto quadrilateral_type = get_cell_type(CellKind::Quadrilateral);
  std::vector<CellType> cell_types = {triangle_type, triangle_type, quadrilateral_type};

  oiseau::mesh::Geometry geometry = oiseau::mesh::Geometry(std::move(x), 2);
  oiseau::mesh::Topology topology = oiseau::mesh::Topology(std::move(conn), std::move(cell_types));
  oiseau::mesh::Mesh mesh(std::move(topology), std::move(geometry));

  auto &g = mesh.geometry();

  std::vector<std::size_t> shape = {g.x().size() / g.dim(), g.dim()};
  auto nodes = xt::adapt(g.x().data(), g.x().size(), xt::no_ownership(), shape);
  auto x_nodes = xt::col(nodes, 0);
  auto y_nodes = xt::col(nodes, 1);

  plt::scoped_interpreter guard{};

  auto [fig, ax] = plt::subplots(1, 1, py::arg("constrained_layout") = true,
                                 "figsize"_a = py::make_tuple(4.0, 2.0));

  ax.scatter(x_nodes, y_nodes, "label"_a = "scatter");
  oiseau::plotting::triplot(ax, mesh);

  plt::show();
}
