// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <span>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

class DGSpace {
 public:
  /// Cells with the same reference shape and polynomial order.
  ///
  /// `nodes(cell, node, coordinate)` contains the physical coordinate of a
  /// cell-local DG node. `cell_indices[cell]` maps a row back to its cell in
  /// the source mesh.
  struct ElementGroup {
    mesh::CellKind kind = mesh::CellKind::Undefined;
    unsigned order = 0;
    std::vector<std::size_t> cell_indices;
    xt::xarray<double> nodes;
  };

  DGSpace(const DGSpace& V) = delete;
  /// Builds a DG space for the mesh's top-dimensional cells.
  ///
  /// `orders` has one positive entry per top-dimensional mesh cell, in mesh
  /// cell order. Lower-dimensional entities, such as boundary edges, are not
  /// assigned an order.
  DGSpace(const mesh::Mesh& mesh, const std::vector<unsigned>& orders);
  DGSpace(DGSpace&& V) = default;
  ~DGSpace() = default;
  DGSpace& operator=(const DGSpace& V) = delete;
  DGSpace& operator=(DGSpace&& V) = delete;

  const mesh::Mesh& mesh() const { return m_mesh; }
  /// Groups indexed by cell shape and polynomial order.
  std::span<const ElementGroup> element_groups() const;
  /// Polynomial orders for top-dimensional cells, in mesh cell order.
  std::span<const unsigned> orders() const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<ElementGroup> m_element_groups;
  const std::vector<unsigned> m_orders;
};

}  // namespace oiseau::dg
