// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_space.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xbuffer_adaptor.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {
namespace {

using GroupKey = std::pair<mesh::CellKind, unsigned>;

struct ActiveCell {
  std::size_t mesh_index = 0;
  mesh::CellType type = nullptr;
  unsigned order = 0;
};

nodal::RefElementType ref_type_for(mesh::CellKind kind) {
  switch (kind) {
  case mesh::CellKind::Interval:
    return nodal::RefElementType::Line;
  case mesh::CellKind::Triangle:
    return nodal::RefElementType::Triangle;
  case mesh::CellKind::Quadrilateral:
    return nodal::RefElementType::Quadrilateral;
  case mesh::CellKind::Tetrahedron:
    return nodal::RefElementType::Tetrahedron;
  case mesh::CellKind::Hexahedron:
    return nodal::RefElementType::Hexahedron;
  default:
    throw std::invalid_argument("DGSpace does not support this mesh cell type");
  }
}

std::vector<ActiveCell> active_cells(std::span<const mesh::CellType> cell_types,
                                     const std::vector<unsigned>& orders, int dimension) {
  const auto count = static_cast<std::size_t>(std::ranges::count_if(
      cell_types, [dimension](const auto* cell) { return cell->dimension() == dimension; }));
  if (orders.size() != count) {
    throw std::invalid_argument(
        "A polynomial order must be provided for every top-dimensional cell");
  }

  std::vector<ActiveCell> cells;
  cells.reserve(count);
  std::size_t order_index = 0;
  for (std::size_t mesh_index = 0; mesh_index < cell_types.size(); ++mesh_index) {
    const auto* type = cell_types[mesh_index];
    if (type->dimension() != dimension) continue;

    const auto order = orders[order_index++];
    if (order == 0) throw std::invalid_argument("Polynomial orders must be greater than zero");
    (void)ref_type_for(type->kind());
    cells.push_back({mesh_index, type, order});
  }
  return cells;
}

std::vector<std::size_t> interpolation_vertices(mesh::CellType type,
                                                std::vector<std::size_t> vertices) {
  if (vertices.size() != static_cast<std::size_t>(type->num_sub_entities(0))) {
    throw std::invalid_argument("Cell connectivity does not match its cell type");
  }
  const auto kind = type->kind();
  if (kind == mesh::CellKind::Quadrilateral || kind == mesh::CellKind::Hexahedron) {
    // Tensor-product reference nodes vary in the first coordinate first;
    // mesh connectivity enumerates tensor-product corners cyclically.
    std::swap(vertices[2], vertices[3]);
    if (kind == mesh::CellKind::Hexahedron) std::swap(vertices[6], vertices[7]);
  }
  return vertices;
}

}  // namespace

DGSpace::DGSpace(const mesh::Mesh& mesh, const std::vector<unsigned>& orders)
    : m_mesh(mesh), m_orders(orders) {
  const auto& topology = mesh.topology();
  const auto& geometry = mesh.geometry();
  const auto cell_types = topology.cell_types();
  const auto connectivity = topology.conn();

  if (cell_types.size() != connectivity.size()) {
    throw std::invalid_argument("Mesh cell types and connectivity counts must match");
  }
  if (geometry.dim() == 0 || geometry.x().size() % geometry.dim() != 0) {
    throw std::invalid_argument("Mesh coordinate data is inconsistent with its dimension");
  }

  const auto dimension = topology.dimension();
  if (dimension < 1 || dimension > 3) {
    throw std::invalid_argument("DGSpace supports meshes with dimensions 1 through 3");
  }
  const auto cells = active_cells(cell_types, orders, dimension);

  std::map<GroupKey, std::size_t> group_sizes;
  for (const auto& cell : cells) ++group_sizes[{cell.type->kind(), cell.order}];

  std::map<GroupKey, std::size_t> group_indices;
  for (const auto& [key, count] : group_sizes) {
    const auto node_count =
        nodal::get_ref_element(ref_type_for(key.first), key.second)->number_of_nodes();
    ElementGroup group{key.first,
                       key.second,
                       {},
                       xt::xarray<double>::from_shape({count, node_count, geometry.dim()})};
    group.cell_indices.reserve(count);
    group_indices.emplace(key, m_element_groups.size());
    m_element_groups.emplace_back(std::move(group));
  }

  const std::array<std::size_t, 2> node_shape = {geometry.x().size() / geometry.dim(),
                                                 geometry.dim()};
  auto mesh_nodes =
      xt::adapt(geometry.x().data(), geometry.x().size(), xt::no_ownership(), node_shape);

  for (const auto& cell : cells) {
    const auto kind = cell.type->kind();
    auto vertices = interpolation_vertices(cell.type, connectivity[cell.mesh_index]);
    if (std::ranges::any_of(vertices,
                            [&node_shape](const auto vertex) { return vertex >= node_shape[0]; })) {
      throw std::out_of_range("Cell references a vertex outside the mesh geometry");
    }

    const auto ref_type = ref_type_for(kind);
    const auto interpolation_element = nodal::get_ref_element(ref_type, 1);
    const auto element = nodal::get_ref_element(ref_type, cell.order);
    const auto vertices_view = xt::view(mesh_nodes, xt::keep(vertices), xt::all());
    const auto interpolation_matrix =
        xt::linalg::dot(interpolation_element->vandermonde(element->r()),
                        xt::linalg::inv(interpolation_element->v()));
    const auto physical_nodes = xt::linalg::dot(interpolation_matrix, vertices_view);

    auto& group = m_element_groups[group_indices.at({kind, cell.order})];
    const auto group_cell_index = group.cell_indices.size();
    xt::view(group.nodes, group_cell_index, xt::all(), xt::all()) = physical_nodes;
    group.cell_indices.push_back(cell.mesh_index);
  }
}

std::span<const DGSpace::ElementGroup> DGSpace::element_groups() const {
  return {m_element_groups};
}

std::span<const unsigned> DGSpace::orders() const { return {m_orders}; }

}  // namespace oiseau::dg
