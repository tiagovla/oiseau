// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_connectivity.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau::dg {
namespace {

std::vector<std::size_t> face_vertices(const mesh::Cell& cell,
                                       const std::vector<std::size_t>& cell_vertices,
                                       std::size_t face_index) {
  const auto local_faces = cell.get_entity_vertices(cell.dimension() - 1);
  if (face_index >= local_faces.size()) throw std::out_of_range("Cell face index is invalid");

  std::vector<std::size_t> vertices;
  vertices.reserve(local_faces[face_index].size());
  for (const auto local_vertex : local_faces[face_index]) {
    if (local_vertex < 0 || static_cast<std::size_t>(local_vertex) >= cell_vertices.size()) {
      throw std::invalid_argument("Cell connectivity does not match its cell type");
    }
    vertices.push_back(cell_vertices[static_cast<std::size_t>(local_vertex)]);
  }
  return vertices;
}

std::vector<std::size_t> vertex_permutation(const std::vector<std::size_t>& local_vertices,
                                            const std::vector<std::size_t>& neighbor_vertices) {
  if (local_vertices.size() != neighbor_vertices.size()) {
    throw std::invalid_argument("Adjacent faces have different numbers of vertices");
  }

  std::vector<std::size_t> permutation;
  permutation.reserve(local_vertices.size());
  for (const auto vertex : local_vertices) {
    const auto position = std::ranges::find(neighbor_vertices, vertex);
    if (position == neighbor_vertices.end()) {
      throw std::invalid_argument("Adjacent faces do not share the same vertices");
    }
    permutation.push_back(static_cast<std::size_t>(position - neighbor_vertices.begin()));
  }
  return permutation;
}

}  // namespace

DGConnectivity::DGConnectivity(mesh::Mesh& mesh) : m_mesh(mesh) {
  auto& topology = mesh.topology();
  topology.calculate_connectivity();

  const auto cell_types = topology.cell_types();
  const auto connectivity = topology.conn();
  const auto cell_neighbors = topology.e_to_e();
  const auto neighbor_faces = topology.e_to_f();
  const auto dimension = topology.dimension();
  if (dimension != 2) {
    throw std::invalid_argument("DGConnectivity currently supports 2D meshes only");
  }

  m_faces.resize(topology.n_cells());
  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    const auto* cell = cell_types[cell_index];
    if (cell->dimension() != dimension) continue;

    const auto local_faces = cell->get_entity_vertices(dimension - 1);
    if (cell_neighbors[cell_index].size() != local_faces.size() ||
        neighbor_faces[cell_index].size() != local_faces.size()) {
      throw std::runtime_error("Mesh face connectivity is inconsistent");
    }
    m_faces[cell_index].resize(local_faces.size());

    for (std::size_t face_index = 0; face_index < local_faces.size(); ++face_index) {
      const auto neighbor_cell = cell_neighbors[cell_index][face_index];
      if (neighbor_cell == cell_index) continue;
      if (neighbor_cell >= cell_types.size()) {
        throw std::out_of_range("Mesh face references an invalid neighboring cell");
      }

      const auto neighbor_face = neighbor_faces[cell_index][face_index];
      const auto* neighbor_type = cell_types[neighbor_cell];
      const auto local_vertices = face_vertices(*cell, connectivity[cell_index], face_index);
      const auto neighbor_vertices =
          face_vertices(*neighbor_type, connectivity[neighbor_cell], neighbor_face);

      auto& connection = m_faces[cell_index][face_index];
      connection.neighbor_cell = neighbor_cell;
      connection.neighbor_face = neighbor_face;
      connection.vertex_permutation = vertex_permutation(local_vertices, neighbor_vertices);
    }
  }
}

std::span<const DGConnectivity::FaceConnection> DGConnectivity::faces(
    std::size_t mesh_cell_index) const {
  if (mesh_cell_index >= m_faces.size()) throw std::out_of_range("Mesh cell index is invalid");
  return m_faces[mesh_cell_index];
}

}  // namespace oiseau::dg
