// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_space.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xbuffer_adaptor.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/element.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

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

  int mesh_dimension = 0;
  std::size_t active_cells = 0;
  for (const auto* cell : cell_types) {
    if (cell == nullptr) throw std::invalid_argument("Mesh cell type cannot be null");
    mesh_dimension = std::max(mesh_dimension, cell->dimension());
    if (cell->dimension() == 2) ++active_cells;
  }
  if (mesh_dimension != 2) {
    throw std::invalid_argument("DGSpace currently supports 2D meshes only");
  }
  if (orders.size() != active_cells) {
    throw std::invalid_argument("A polynomial order must be provided for every 2D mesh cell");
  }

  m_elements.reserve(active_cells);
  m_cell_indices.reserve(active_cells);

  std::array<std::size_t, 2> shape = {geometry.x().size() / geometry.dim(), geometry.dim()};
  auto nodes = xt::adapt(geometry.x().data(), geometry.x().size(), xt::no_ownership(), shape);

  std::size_t order_index = 0;
  for (std::size_t i = 0; i < cell_types.size(); ++i) {
    const auto& cell_type = cell_types[i];
    if (cell_type->dimension() != 2) continue;
    mesh::CellKind kind = cell_type->kind();

    nodal::RefElementType ref_type;
    switch (kind) {
    case mesh::CellKind::Triangle:
      ref_type = nodal::RefElementType::Triangle;
      break;
    case mesh::CellKind::Quadrilateral:
      ref_type = nodal::RefElementType::Quadrilateral;
      break;
    default:
      throw std::runtime_error("Unsupported cell type");
    }

    const auto& cell_vertices = connectivity[i];
    const std::size_t expected_vertices = kind == mesh::CellKind::Triangle ? 3 : 4;
    if (cell_vertices.size() != expected_vertices) {
      throw std::invalid_argument("Cell connectivity does not match its cell type");
    }
    if (std::any_of(cell_vertices.begin(), cell_vertices.end(),
                    [&shape](const auto vertex) { return vertex >= shape[0]; })) {
      throw std::out_of_range("Cell references a vertex outside the mesh geometry");
    }
    if (orders[order_index] == 0) {
      throw std::invalid_argument("Polynomial orders must be greater than zero");
    }

    auto interpolation_vertices = cell_vertices;
    if (kind == mesh::CellKind::Quadrilateral) {
      // The reference nodes are tensor ordered (bottom-left, bottom-right,
      // top-left, top-right), while the mesh cell uses cyclic order.
      std::swap(interpolation_vertices[2], interpolation_vertices[3]);
    }

    auto interp_elem = nodal::get_ref_element(ref_type, 1);
    auto ref_elem = nodal::get_ref_element(ref_type, orders[order_index]);
    auto x_view = xt::view(nodes, xt::keep(interpolation_vertices), xt::all());

    auto inv_v = xt::linalg::inv(interp_elem->v());
    auto v = interp_elem->vandermonde(ref_elem->r());
    auto interp_x = xt::linalg::dot(xt::linalg::dot(v, inv_v), x_view);

    m_elements.emplace_back(ref_elem, interp_x);
    m_cell_indices.emplace_back(i);
    ++order_index;
  }
  // TODO(tiagovla): clean up this mess, introduce proper api
}
std::span<const nodal::Element> DGSpace::elements() const { return {m_elements}; }
std::span<const unsigned> DGSpace::orders() const { return {m_orders}; }
std::span<const std::size_t> DGSpace::cell_indices() const { return {m_cell_indices}; }

}  // namespace oiseau::dg
