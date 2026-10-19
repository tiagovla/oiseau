// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_space.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/plotting/triplot.hpp"

namespace {

const char* cell_name(oiseau::mesh::CellKind kind) {
  using oiseau::mesh::CellKind;
  switch (kind) {
  case CellKind::Interval:
    return "line";
  case CellKind::Triangle:
    return "triangle";
  case CellKind::Quadrilateral:
    return "quadrilateral";
  case CellKind::Tetrahedron:
    return "tetrahedron";
  case CellKind::Hexahedron:
    return "hexahedron";
  default:
    return "unknown";
  }
}

int topological_dimension(std::span<const oiseau::mesh::CellType> cell_types) {
  int dimension = 0;
  for (const auto* cell : cell_types) dimension = std::max(dimension, cell->dimension());
  return dimension;
}

void print_group(std::size_t group_index, const oiseau::dg::DGSpace::ElementGroup& group) {
  const auto& shape = group.nodes.shape();
  if (shape.size() != 3) {
    throw std::runtime_error("DG group coordinates must have shape [cells, nodes, dimensions]");
  }

  const auto cell_count = shape[0];
  const auto node_count = shape[1];
  const auto coordinate_count = shape[2];
  std::cout << "Group " << group_index << ": " << cell_name(group.kind) << ", order " << group.order
            << "\n"
            << "  nodes: [" << cell_count << " cells x " << node_count << " nodes x "
            << coordinate_count << " coordinates]\n";

  // A group row is a cell. cell_indices maps that row to the original mesh.
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    std::cout << "  group cell " << cell << " (mesh cell " << group.cell_indices[cell] << "):\n";
    for (std::size_t node = 0; node < node_count; ++node) {
      std::cout << "    node " << node << " = (";
      for (std::size_t coordinate = 0; coordinate < coordinate_count; ++coordinate) {
        if (coordinate != 0) std::cout << ", ";
        std::cout << group.nodes(cell, node, coordinate);
      }
      std::cout << ")\n";
    }
  }
}

std::vector<std::vector<double>> collect_node_coordinates(
    std::span<const oiseau::dg::DGSpace::ElementGroup> groups, std::size_t dimension) {
  std::vector<std::vector<double>> coordinates(dimension);
  for (const auto& group : groups) {
    const auto& shape = group.nodes.shape();
    for (std::size_t cell = 0; cell < shape[0]; ++cell) {
      for (std::size_t node = 0; node < shape[1]; ++node) {
        for (std::size_t coordinate = 0; coordinate < shape[2]; ++coordinate) {
          coordinates[coordinate].push_back(group.nodes(cell, node, coordinate));
        }
      }
    }
  }
  return coordinates;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace oiseau;

  const std::filesystem::path mesh_path =
      argc > 1 ? argv[1] : "demo/meshes/mesh2d_two_triangles_one_quad.msh";
  auto mesh = io::gmsh_read_from_path(mesh_path);
  const auto cell_types = mesh.topology().cell_types();

  const auto mesh_dimension = topological_dimension(cell_types);

  // One order for each top-dimensional cell: triangle, triangle, then quad.
  const std::vector<unsigned> orders{1, 1, 3};

  // DGSpace groups physical nodal coordinates by cell kind and polynomial
  // order. Each group stores nodes as [cells, nodes, dimensions].
  const dg::DGSpace space(mesh, orders);
  const auto groups = space.element_groups();

  std::cout << "DGSpace demo\n"
            << "Mesh: " << mesh_path << "\n"
            << "Top-dimensional cells: " << orders.size() << " (dimension " << mesh_dimension
            << ")\n"
            << "Lower-dimensional entities skipped: " << cell_types.size() - orders.size()
            << "\nGroups: " << groups.size() << "\n\n";

  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index)
    print_group(group_index, groups[group_index]);

  // DG nodes are local to each cell, so shared physical points may appear more
  // than once across groups or neighboring cells.
  const auto node_coordinates = collect_node_coordinates(groups, mesh.geometry().dim());

  // triplot is a two-dimensional view. Coordinate generation and reporting
  // above work for line, surface, and volume elements.
  if (mesh_dimension == 2) {
    plt::scoped_interpreter interpreter{};
    auto [figure, axes] =
        plt::subplots(1, 1, "figsize"_a = py::make_tuple(8.0, 6.0), "constrained_layout"_a = true);
    plotting::triplot(axes, mesh);
    axes.scatter(node_coordinates[0], node_coordinates[1], "s"_a = 9, "color"_a = "tab:red",
                 "label"_a = "DG nodes");
    axes.set_aspect("equal");
    axes.set_xlabel("x");
    axes.set_ylabel("y");
    axes.set_title("DGSpace nodes");
    axes.legend();
    plt::show();
  }

  return 0;
}
