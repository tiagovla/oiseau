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

/**
 * @brief Stores geometric factors for top-dimensional cells in a mesh.
 * @details The mesh is referenced, not copied, and must outlive this object.
 * Geometry entries are indexed by source mesh cell. Lower-dimensional cells
 * have no entry.
 */
class Geometry {
 public:
  /** @brief Geometric data for one cell face. */
  struct FaceGeometry {
    /** @brief Outward unit normal, padded to three coordinates. */
    std::array<double, 3> outward_normal{};
    /**
     * @brief Physical face measure.
     * @details This is face length in 2D and face area in 3D. In 1D, it is 1.
     */
    double measure = 0.0;
    /**
     * @brief Ratio of physical to reference face measure.
     * @details Multiply reference face weights by this value to scale them.
     */
    double surface_jacobian = 0.0;
  };

  /** @brief Geometric data for one top-dimensional mesh cell. */
  struct CellGeometry {
    /** @brief Index of this cell in the source mesh. */
    std::size_t mesh_cell_index = 0;
    /** @brief Topological dimension of this cell. */
    std::size_t dimension = 0;
    /**
     * @brief Jacobian mapping reference coordinates to physical coordinates.
     * @details Padded row-major 3 by 3 storage. Entry (row, column) is at
     * row * 3 + column.
     */
    std::array<double, 9> jacobian{};
    /**
     * @brief Inverse Jacobian in padded row-major 3 by 3 storage.
     * @details Entry (row, column) is at row * 3 + column.
     */
    std::array<double, 9> inverse_jacobian{};
    /**
     * @brief Signed Jacobian determinant.
     * @details Use its absolute value to scale volume integrals.
     */
    double determinant = 0.0;
    /** @brief Face geometry in the cell type's local face order. */
    std::vector<FaceGeometry> faces;
  };

  /**
   * @brief Computes geometric factors for the mesh's top-dimensional cells.
   * @param mesh Source mesh. It must outlive this object.
   * @throws std::invalid_argument If dimensions, coordinates, or cell geometry
   * are invalid, or if a cell mapping is non-affine.
   * @throws std::out_of_range If a cell references a missing mesh vertex.
   */
  explicit Geometry(const mesh::Mesh& mesh);

  /** @brief Returns the source mesh. */
  const mesh::Mesh& mesh() const { return m_mesh; }
  /**
   * @brief Checks whether geometric data exists for a mesh cell.
   * @param mesh_cell_index Source mesh cell index.
   * @return True if the cell is top-dimensional and has geometry data.
   */
  bool has_cell(std::size_t mesh_cell_index) const;
  /**
   * @brief Returns geometric data for a mesh cell.
   * @param mesh_cell_index Source mesh cell index.
   * @return The cell's geometric data.
   * @throws std::out_of_range If the index is invalid or the cell has no data.
   */
  const CellGeometry& cell(std::size_t mesh_cell_index) const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<std::optional<CellGeometry>> m_cells;
};

}  // namespace oiseau::dg
