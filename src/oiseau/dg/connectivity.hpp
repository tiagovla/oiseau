// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

/**
 * @brief Stores face adjacency for the top-dimensional cells of a 1D, 2D, or 3D mesh.
 * @details Mesh topology supplies the neighbor indices. Interior faces also
 * store a vertex permutation so a numerical flux can align data from both
 * sides of the interface. Connections follow the cell type's local face order
 * and are indexed by the original mesh cell index.
 */
class DGConnectivity {
 public:
  /** @brief Sentinel for missing neighbor cell and face indices. */
  static constexpr std::size_t no_neighbor = std::numeric_limits<std::size_t>::max();

  /** @brief Connection from one local face to its neighbor. */
  struct FaceConnection {
    /** @brief Mesh index of the neighboring cell, or @c no_neighbor on a boundary. */
    std::size_t neighbor_cell = no_neighbor;
    /** @brief Local face index on the neighboring cell, or @c no_neighbor on a boundary. */
    std::size_t neighbor_face = no_neighbor;
    /**
     * @brief Maps local face vertex positions to positions on the neighboring face.
     * @details @c vertex_permutation[i] gives the position of this face's vertex
     * @c i in the neighbor's face ordering. A boundary face has an empty permutation.
     */
    std::vector<std::size_t> vertex_permutation;

    /** @return True if this face has no neighboring cell. */
    bool is_boundary() const { return neighbor_cell == no_neighbor; }
  };

  /**
   * @brief Calculates mesh adjacency and builds DG face connections.
   * @details Updates the source topology through
   * @c Topology::calculate_connectivity(). Rebuild after changing the mesh topology.
   * @param mesh Source mesh. It must outlive this object.
   * @throws std::invalid_argument If the mesh dimension is outside 1 through 3 or its connectivity
   * is invalid.
   * @throws std::out_of_range If adjacency references an invalid neighbor cell or face.
   * @throws std::runtime_error If the topology's adjacency counts differ from the cell's face
   * count.
   */
  explicit DGConnectivity(mesh::Mesh& mesh);

  /** @return A const reference to the source mesh. */
  const mesh::Mesh& mesh() const { return m_mesh; }

  /**
   * @brief Returns the face connections for a mesh cell in local face order.
   * @param mesh_cell_index Index in the source mesh.
   * @return A view of the stored connections. Lower-dimensional cells return an empty span.
   * @throws std::out_of_range If the mesh cell index is invalid.
   */
  std::span<const FaceConnection> faces(std::size_t mesh_cell_index) const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<std::vector<FaceConnection>> m_faces;
};

}  // namespace oiseau::dg
