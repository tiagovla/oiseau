// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

/// Constant geometric factors for affine cells of every supported reference type.
///
/// The Jacobian maps reference coordinates to physical coordinates. Matrices
/// use a 3-column padded row-major layout: entry (row, column) is at
/// `row * 3 + column`. Coordinates must lie in the subspace spanned by the
/// first `dimension()` geometry components.
class DGGeometry {
 public:
  struct FaceGeometry {
    std::array<double, 3> outward_normal{};
    double measure = 0.0;
    /// Physical face measure divided by the corresponding reference measure.
    double surface_jacobian = 0.0;
  };

  struct CellGeometry {
    std::size_t mesh_cell_index = 0;
    std::size_t dimension = 0;
    /// Row-major `dimension x dimension` matrix d(x) / d(r).
    std::array<double, 9> jacobian{};
    /// Row-major inverse of `jacobian`.
    std::array<double, 9> inverse_jacobian{};
    /// Signed determinant; use its absolute value for volume integration.
    double determinant = 0.0;
    std::vector<FaceGeometry> faces;
  };

  /// Builds geometry for all top-dimensional affine mesh cells.
  ///
  /// Lines, triangles, quadrilaterals, tetrahedra, and hexahedra are
  /// supported. Tensor-product cells must be affine: quadrilaterals must be
  /// parallelograms and hexahedra must be parallelepipeds.
  explicit DGGeometry(const mesh::Mesh& mesh);

  const mesh::Mesh& mesh() const { return m_mesh; }
  bool has_cell(std::size_t mesh_cell_index) const;
  const CellGeometry& cell(std::size_t mesh_cell_index) const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<std::optional<CellGeometry>> m_cells;
};

}  // namespace oiseau::dg
