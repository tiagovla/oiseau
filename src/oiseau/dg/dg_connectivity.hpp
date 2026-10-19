// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
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

/// Face adjacency for the top-dimensional cells of a two-dimensional mesh.
///
/// The mesh topology supplies the neighboring cell and face. DGConnectivity
/// additionally records how the vertices on an interior face correspond, so
/// a numerical flux can align data from both sides of the interface.
class DGConnectivity {
 public:
  static constexpr std::size_t no_neighbor = std::numeric_limits<std::size_t>::max();

  struct FaceConnection {
    std::size_t neighbor_cell = no_neighbor;
    std::size_t neighbor_face = no_neighbor;
    /// `vertex_permutation[i]` is the position of this face's vertex `i` in
    /// the neighboring face's vertex ordering.
    std::vector<std::size_t> vertex_permutation;

    bool is_boundary() const { return neighbor_cell == no_neighbor; }
  };

  /// Calculates mesh adjacency and derives the DG face orientation data.
  ///
  /// The supplied mesh is updated through `Topology::calculate_connectivity()`.
  /// This first implementation supports top-dimensional 2D cells.
  explicit DGConnectivity(mesh::Mesh& mesh);

  const mesh::Mesh& mesh() const { return m_mesh; }
  /// Connections for local faces of a mesh cell. Lower-dimensional mesh
  /// entities have no entries.
  std::span<const FaceConnection> faces(std::size_t mesh_cell_index) const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<std::vector<FaceConnection>> m_faces;
};

}  // namespace oiseau::dg
