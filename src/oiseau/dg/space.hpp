// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <span>
#include <vector>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

/**
 * @brief Stores discontinuous Galerkin nodes grouped by cell kind and order.
 * @details The mesh is referenced, not copied, and must outlive this space.
 */
class DGSpace {
 public:
  /** @brief Stores cells that share a cell kind and polynomial order. */
  struct ElementGroup {
    /** @brief Reference-cell kind shared by every cell in the group. */
    mesh::CellKind kind = mesh::CellKind::Undefined;
    /** @brief Polynomial order shared by every cell in the group. */
    unsigned order = 0;
    /** @brief Source mesh indices for cells in this group. */
    std::vector<std::size_t> cell_indices;
    /**
     * @brief Physical coordinates of the group nodes.
     * @details Shape is [cells, nodes per cell, coordinate dimensions].
     */
    xt::xarray<double> nodes;
  };

  DGSpace(const DGSpace& V) = delete;
  /**
   * @brief Builds groups using one order per top-dimensional mesh cell.
   * @param mesh Source mesh. It must outlive this space.
   * @param orders Positive polynomial orders in top-dimensional cell order.
   * Lower-dimensional cells, such as boundary edges, have no order.
   * @throws std::invalid_argument If the order count, mesh dimensions,
   * geometry, or cell kinds are invalid.
   * @throws std::out_of_range If a cell references a missing mesh vertex.
   */
  DGSpace(const mesh::Mesh& mesh, const std::vector<unsigned>& orders);
  /**
   * @brief Builds groups using one order for every top-dimensional cell.
   * @param mesh Source mesh. It must outlive this space.
   * @param order Positive polynomial order assigned to each top-dimensional cell.
   * @throws std::invalid_argument If the order is zero or the mesh is invalid.
   * @throws std::out_of_range If a cell references a missing mesh vertex.
   */
  DGSpace(const mesh::Mesh& mesh, unsigned order);
  DGSpace(DGSpace&& V) = default;
  ~DGSpace() = default;
  DGSpace& operator=(const DGSpace& V) = delete;
  DGSpace& operator=(DGSpace&& V) = delete;

  /** @brief Returns the source mesh. */
  const mesh::Mesh& mesh() const { return m_mesh; }
  /** @brief Returns all element groups. */
  std::span<const ElementGroup> element_groups() const;
  /**
   * @brief Finds a group by cell kind and polynomial order.
   * @param kind Cell kind to find.
   * @param order Polynomial order to find.
   * @return Matching group, or nullptr if no group exists.
   */
  const ElementGroup* element_group(mesh::CellKind kind, unsigned order) const;
  /** @brief Returns orders for top-dimensional cells in mesh order. */
  std::span<const unsigned> orders() const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<ElementGroup> m_element_groups;
  const std::vector<unsigned> m_orders;
};

}  // namespace oiseau::dg
