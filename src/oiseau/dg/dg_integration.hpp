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
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/dg_space.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

/**
 * @brief Precomputed volume and face integration operators for affine 1D, 2D, and 3D DG cells.
 * @details The rules follow the supplied DGSpace cell groups and can be reused by
 * different equation systems on the same space.
 */
class DGIntegration {
 public:
  /** @brief Array type for integration operators and quadrature data. */
  using Matrix = xt::xarray<double>;

  /** @brief Reference-cell volume operators shared by one cell shape and polynomial order. */
  struct ReferenceRule {
    /** @brief Reference mass matrix. */
    Matrix mass;
    /** @brief Inverse reference mass matrix. */
    Matrix inverse_mass;
    /** @brief Reference weak stiffness matrix in the @c r direction. */
    Matrix stiffness_r;
    /** @brief Reference weak stiffness matrix in the @c s direction. */
    Matrix stiffness_s;
    /** @brief Reference weak stiffness matrix in the @c t direction. */
    Matrix stiffness_t;
    /** @brief Inverse mass times reference stiffness in the @c r direction. */
    Matrix weak_derivative_r;
    /** @brief Inverse mass times reference stiffness in the @c s direction. */
    Matrix weak_derivative_s;
    /** @brief Inverse mass times reference stiffness in the @c t direction. */
    Matrix weak_derivative_t;
    /** @brief Strong nodal derivative matrix in the @c r direction. */
    Matrix derivative_r;
    /** @brief Strong nodal derivative matrix in the @c s direction. */
    Matrix derivative_s;
    /** @brief Strong nodal derivative matrix in the @c t direction. */
    Matrix derivative_t;
  };

  /**
   * @brief Weights and nodes on the canonical point, line, triangle, or square.
   * @details A point has a single unit weight and no local coordinates. For other shapes,
   * the rule integrates degree @f$2p+1@f$ exactly, per coordinate on the square.
   * Here @f$p@f$ is the higher of the adjacent cell orders on an interior face and
   * the local cell order on a boundary face.
   */
  struct QuadratureRule {
    /** @brief Coordinates of the quadrature points on the canonical face. */
    Matrix nodes;
    /** @brief Weights in quadrature point order. */
    Matrix weights;
  };

  /** @brief Face operators shared by cell shape, order, local face, quadrature, and orientation. */
  struct ReferenceFaceRule {
    /** @brief Canonical face quadrature used by the operators. */
    std::shared_ptr<const QuadratureRule> quadrature;
    /** @brief Trace interpolation matrix with shape [quadrature point, cell node]. */
    Matrix trace;
    /**
     * @brief Reference lift @f$M^{-1}T^TW@f$, stored as [quadrature point, test node].
     * @details @f$W@f$ uses canonical face weights. Physical scaling belongs to FaceRule.
     */
    Matrix lift;
  };

  /** @brief Integration data for one physical cell face. */
  struct FaceRule {
    /** @brief The x and y components of the outward unit normal. */
    std::array<double, 2> outward_normal{};
    /** @brief The z component of the outward unit normal. */
    double normal_z = 0.0;
    /** @brief Mesh index of the neighbor cell, or DGConnectivity::no_neighbor on a boundary. */
    std::size_t neighbor_cell = DGConnectivity::no_neighbor;
    /** @brief Reference face rule for the local cell. */
    std::shared_ptr<const ReferenceFaceRule> minus;
    /**
     * @brief Reference face rule for the neighbor cell, or @c nullptr on a boundary.
     * @details On interior faces, both traces evaluate at the same physical points
     * in the same order. The plus trace uses the neighbor's nodes.
     */
    std::shared_ptr<const ReferenceFaceRule> plus;
    /** @brief Ratio of physical face measure to canonical face measure. */
    double weight_scale = 0.0;
    /** @brief The face weight scale divided by @f$|\det J|@f$ to scale the reference lift. */
    double lift_scale = 0.0;

    /** @return Number of quadrature points on the face. */
    std::size_t point_count() const { return minus->quadrature->weights.size(); }
    /**
     * @brief Returns a quadrature weight scaled by the physical face measure ratio.
     * @param q Quadrature point index.
     */
    double physical_weight(std::size_t q) const {
      return weight_scale * minus->quadrature->weights(q);
    }
    /**
     * @brief Returns an entry of the physical lift for the local cell.
     * @param q Quadrature point index.
     * @param node Local test node index.
     */
    double physical_lift(std::size_t q, std::size_t node) const {
      return lift_scale * minus->lift(q, node);
    }

    /** @return True if the face has no neighbor cell. */
    bool is_boundary() const { return neighbor_cell == DGConnectivity::no_neighbor; }
  };

  /** @brief Volume and face integration data for one mesh cell. */
  struct CellRule {
    /** @brief Index of the DGSpace group, or DGConnectivity::no_neighbor for an inactive cell. */
    std::size_t group_index = DGConnectivity::no_neighbor;
    /** @brief Index of this cell within its DGSpace group. */
    std::size_t group_cell_index = 0;
    /** @brief Shared reference volume operators, or @c nullptr for an inactive cell. */
    std::shared_ptr<const ReferenceRule> reference;
    /** @brief Absolute Jacobian determinant used to scale volume integrals. */
    double volume_jacobian = 0.0;
    /** @brief Upper-left 2 by 2 block of the inverse Jacobian in row-major order. */
    std::array<double, 4> inverse_jacobian{};
    /** @brief Padded row-major 3 by 3 inverse Jacobian @f$\partial(r,s,t)/\partial(x,y,z)@f$. */
    std::array<double, 9> spatial_inverse_jacobian{};

    /** @return Number of cell nodes, or zero for an inactive cell. */
    std::size_t node_count() const { return reference ? reference->mass.shape()[0] : 0; }
    /** @return Physical mass matrix, or an empty matrix for an inactive cell. */
    Matrix mass() const;
    /** @return Inverse physical mass matrix, or an empty matrix for an inactive cell. */
    Matrix inverse_mass() const;
    /**
     * @brief Returns the physical weak stiffness matrix in the x direction.
     * @details @f$S_x(i,j) = \int (\partial_x\phi_i)\phi_j = (D_x^TM)(i,j)@f$.
     * Rows are test nodes and columns are trial nodes. The y and z matrices follow
     * the same convention.
     * @return The stiffness matrix, or an empty matrix for an inactive cell.
     */
    Matrix stiffness_x() const;
    /** @return Physical y stiffness matrix, or an empty matrix if y is absent. */
    Matrix stiffness_y() const;
    /** @return Physical z stiffness matrix, or an empty matrix if z is absent. */
    Matrix stiffness_z() const;
    /** @brief Physical weak derivative @f$M^{-1}S_x@f$, stored as [test node, trial node]. */
    std::vector<double> weak_derivative_x;
    /** @brief Physical weak derivative @f$M^{-1}S_y@f$, stored as [test node, trial node]. */
    std::vector<double> weak_derivative_y;
    /** @brief Physical weak derivative @f$M^{-1}S_z@f$, stored as [test node, trial node]. */
    std::vector<double> weak_derivative_z;
    /** @brief Face rules in local face order. */
    std::vector<FaceRule> faces;
  };

  /**
   * @brief Builds integration rules for the cells in a DG space.
   * @param mesh Source mesh.
   * @param space DG cell groups and polynomial orders.
   * @param connectivity Face adjacency and orientation for the source mesh.
   * @param geometry Cell and face geometry for the source mesh.
   * @throws std::invalid_argument If the inputs refer to different meshes, a neighboring cell
   * is missing from the space, or a polynomial order is too large for face quadrature.
   */
  DGIntegration(const mesh::Mesh& mesh, const DGSpace& space, const DGConnectivity& connectivity,
                const DGGeometry& geometry);

  /**
   * @return Rules in mesh cell order. Cells outside the DGSpace have empty matrices.
   */
  std::span<const CellRule> cells() const { return m_cells; }
  /**
   * @brief Returns the integration rule for a mesh cell.
   * @param mesh_cell_index Index in the source mesh.
   * @return The cell rule, with empty matrices for a cell outside the DGSpace.
   * @throws std::out_of_range If the mesh cell index is invalid.
   */
  const CellRule& cell(std::size_t mesh_cell_index) const;

 private:
  std::vector<CellRule> m_cells;
};

}  // namespace oiseau::dg
