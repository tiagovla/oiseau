// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/dg_space.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

/// Precomputed volume and face integration operators for affine 2D DG cells.
///
/// The initial implementation supports triangles and parallelograms. It owns
/// no mesh data; its rules follow the supplied DGSpace cell groups and can be
/// reused by different equation systems on the same space.
class DGIntegration {
 public:
  using Matrix = xt::xarray<double>;

  /// Reference-cell volume operators shared by one (cell shape, order) group.
  struct ReferenceRule {
    Matrix mass;
    Matrix inverse_mass;
    Matrix stiffness_r;
    Matrix stiffness_s;
    Matrix weak_derivative_r;
    Matrix weak_derivative_s;
  };

  struct FaceRule {
    std::array<double, 2> outward_normal{};
    std::size_t neighbor_cell = DGConnectivity::no_neighbor;
    /// Traces have shape [quadrature point, cell node]. On interior faces,
    /// both traces evaluate at the same physical points in the same order;
    /// plus_trace uses the neighbor's nodes and is empty on boundaries.
    Matrix minus_trace;
    Matrix plus_trace;
    /// Positive Gauss weights including physical edge length scaling.
    /// The common rule integrates degree 2 * max(p_minus, p_plus) + 1 exactly.
    std::vector<double> physical_weights;
    /// M^-1 T_minus^T W, stored as [quadrature point, test node].
    std::vector<double> lift;

    bool is_boundary() const { return neighbor_cell == DGConnectivity::no_neighbor; }
  };

  struct CellRule {
    std::size_t group_index = DGConnectivity::no_neighbor;
    std::size_t group_cell_index = 0;
    std::shared_ptr<const ReferenceRule> reference;
    double volume_jacobian = 0.0;
    /// Row-major 2 x 2 matrix d(r,s) / d(x,y).
    std::array<double, 4> inverse_jacobian{};

    std::size_t node_count() const { return reference ? reference->mass.shape()[0] : 0; }
    /// Materialize physical matrices on demand; inactive cells return empty matrices.
    Matrix mass() const;
    Matrix inverse_mass() const;
    /// Weak stiffness: S_x(i,j) = integral (d_x phi_i) phi_j = (D_x^T M)(i,j).
    /// Rows are test nodes, columns are trial nodes; similarly for S_y.
    Matrix stiffness_x() const;
    Matrix stiffness_y() const;
    /// M^-1 S_x and M^-1 S_y, stored as [test node, trial node].
    std::vector<double> weak_derivative_x;
    std::vector<double> weak_derivative_y;
    std::vector<FaceRule> faces;
  };

  DGIntegration(const mesh::Mesh& mesh, const DGSpace& space, const DGConnectivity& connectivity,
                const DGGeometry& geometry);

  /// Rules in mesh cell order. Cells outside the DGSpace have empty matrices.
  std::span<const CellRule> cells() const { return m_cells; }
  const CellRule& cell(std::size_t mesh_cell_index) const;

 private:
  std::vector<CellRule> m_cells;
};

}  // namespace oiseau::dg
