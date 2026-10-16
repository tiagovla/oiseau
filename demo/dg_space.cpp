// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_space.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "oiseau/dg/nodal/element.hpp"
#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/plotting/triplot.hpp"

int main(int argc, char** argv) {
  using namespace oiseau;

  const std::filesystem::path mesh_path =
      argc > 1 ? argv[1] : "demo/meshes/mesh2d_split_triangles_quads.msh";
  auto mesh = io::gmsh_read_from_path(mesh_path);
  const auto cell_types = mesh.topology().cell_types();

  // DGSpace expects one polynomial order for each 2D cell. Use a higher order
  // on triangles and a lower order on quadrilaterals; boundary entities do not
  // get entries in this vector because DGSpace skips them.
  std::vector<unsigned> orders;
  for (const auto* type : cell_types) {
    if (type->dimension() != 2) continue;
    switch (type->kind()) {
    case mesh::CellKind::Triangle:
      orders.push_back(4);
      break;
    case mesh::CellKind::Quadrilateral:
      orders.push_back(2);
      break;
    default:
      throw std::runtime_error("DGSpace demo supports triangles and quadrilaterals only");
    }
  }

  // Constructing DGSpace creates one nodal element per 2D mesh cell at the
  // corresponding order. The physical nodal coordinates are interpolated from
  // each cell's mesh vertices.
  const dg::DGSpace space(mesh, orders);
  const auto elements = space.elements();
  const auto source_cells = space.cell_indices();

  std::size_t triangles = 0;
  std::size_t quadrilaterals = 0;
  for (const auto* type : cell_types) {
    if (type->dimension() != 2) continue;
    if (type->kind() == mesh::CellKind::Triangle) ++triangles;
    if (type->kind() == mesh::CellKind::Quadrilateral) ++quadrilaterals;
  }

  std::cout << "DGSpace demo\n"
            << "Mesh: " << mesh_path << "\n"
            << "2D cells: " << elements.size() << " (" << triangles << " triangles, "
            << quadrilaterals << " quadrilaterals)\n"
            << "Lower-dimensional entities skipped: " << cell_types.size() - elements.size()
            << "\n\n";

  // cell_indices() is the mapping from DG element number back to the original
  // mesh cell number. The matching order is available from orders().
  const std::size_t examples = std::min<std::size_t>(elements.size(), 5);
  for (std::size_t element_index = 0; element_index < examples; ++element_index) {
    const auto& element = elements[element_index];
    const auto& nodes = element.nodes();
    const auto& shape = nodes.shape();
    if (shape.size() != 2 || shape[1] < 2) {
      throw std::runtime_error("DG element coordinates have an unexpected shape");
    }

    std::cout << "DG element " << element_index << " <- mesh cell " << source_cells[element_index]
              << ": order " << element.order() << ", " << element.reference().number_of_nodes()
              << " reference nodes, physical node array [" << shape[0] << " x " << shape[1]
              << "]\n";
    for (std::size_t node = 0; node < shape[0]; ++node) {
      std::cout << "  node " << node << " at (" << nodes(node, 0) << ", " << nodes(node, 1)
                << ")\n";
    }
  }
  if (elements.size() > examples) {
    std::cout << "... " << elements.size() - examples << " more DG elements\n";
  }

  // DG nodes are local to each element. Adjacent elements can therefore have
  // separate degrees of freedom at the same physical coordinate; plotting all
  // of them makes that cell-local nodal layout visible on the mesh.
  std::vector<double> node_x;
  std::vector<double> node_y;
  for (const auto& element : elements) {
    const auto& nodes = element.nodes();
    for (std::size_t node = 0; node < nodes.shape()[0]; ++node) {
      node_x.push_back(nodes(node, 0));
      node_y.push_back(nodes(node, 1));
    }
  }

  // This visualization is only of the DGSpace mesh and its nodal points.
  plt::scoped_interpreter interpreter{};
  auto [figure, axes] =
      plt::subplots(1, 1, "figsize"_a = py::make_tuple(8.0, 6.0), "constrained_layout"_a = true);
  plotting::triplot(axes, mesh);
  axes.scatter(node_x, node_y, "s"_a = 9, "color"_a = "tab:red", "label"_a = "DG nodes");
  axes.set_aspect("equal");
  axes.set_xlabel("x");
  axes.set_ylabel("y");
  axes.set_title("DGSpace: order-4 triangle and order-2 quadrilateral nodes");
  axes.legend();
  plt::show();

  return 0;
}
