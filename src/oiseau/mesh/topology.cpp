// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/topology.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "oiseau/mesh/cell.hpp"

using namespace oiseau::mesh;

Topology::Topology() = default;
Topology::~Topology() = default;

Topology::Topology(std::vector<std::vector<std::size_t>>&& conn, std::vector<CellType>&& cell_types)
    : m_conn(std::move(conn)), m_cell_types(std::move(cell_types)) {};

std::span<CellType> Topology::cell_types() { return m_cell_types; };
std::span<const CellType> Topology::cell_types() const { return m_cell_types; };

std::span<std::vector<std::size_t>> Topology::conn() { return m_conn; };
std::span<const std::vector<std::size_t>> Topology::conn() const { return m_conn; };
std::span<std::vector<std::size_t>> Topology::e_to_e() { return m_e_to_e; };
std::span<std::vector<std::size_t>> Topology::e_to_f() { return m_e_to_f; };

std::size_t Topology::n_cells() const { return m_conn.size(); }

void Topology::calculate_connectivity() {
  if (m_conn.size() != m_cell_types.size()) {
    throw std::invalid_argument("Cell connectivity and cell type counts must match");
  }

  int mesh_dimension = 0;
  for (const auto* cell : m_cell_types) {
    if (cell == nullptr) throw std::invalid_argument("Cell type cannot be null");
    mesh_dimension = std::max(mesh_dimension, cell->dimension());
  }
  if (mesh_dimension > 2) {
    throw std::invalid_argument("Face connectivity currently supports 2D meshes only");
  }

  // Each canonical face key records the first cell/local-face that uses it.
  // A second use creates the two-way adjacency; a third is non-manifold.
  std::map<std::vector<std::size_t>, std::pair<std::size_t, std::size_t>> face_owner;
  std::vector<std::vector<std::size_t>> e_to_e(m_conn.size());
  std::vector<std::vector<std::size_t>> e_to_f(m_conn.size());

  for (std::size_t cell_index = 0; cell_index < m_conn.size(); ++cell_index) {
    const auto* cell = m_cell_types[cell_index];
    // Gmsh meshes commonly include lower-dimensional boundary entities. They
    // are not volume cells and therefore do not participate in this adjacency.
    if (cell->dimension() != mesh_dimension) continue;

    const auto face_vertices = cell->get_entity_vertices(1);
    const auto& cell_vertices = m_conn[cell_index];
    e_to_e[cell_index].resize(face_vertices.size(), cell_index);
    e_to_f[cell_index].resize(face_vertices.size());

    for (std::size_t face_index = 0; face_index < face_vertices.size(); ++face_index) {
      std::vector<std::size_t> key;
      key.reserve(face_vertices[face_index].size());
      for (const int local_vertex : face_vertices[face_index]) {
        if (local_vertex < 0 || static_cast<std::size_t>(local_vertex) >= cell_vertices.size()) {
          throw std::invalid_argument("Cell connectivity does not match its cell type");
        }
        key.push_back(cell_vertices[static_cast<std::size_t>(local_vertex)]);
      }
      std::ranges::sort(key);
      e_to_f[cell_index][face_index] = face_index;

      const auto [owner, inserted] =
          face_owner.emplace(std::move(key), std::make_pair(cell_index, face_index));
      if (!inserted) {
        const auto [other_cell, other_face] = owner->second;
        if (other_cell == cell_index || e_to_e[other_cell][other_face] != other_cell) {
          throw std::invalid_argument("Non-manifold face is shared by more than two cells");
        }
        e_to_e[cell_index][face_index] = other_cell;
        e_to_f[cell_index][face_index] = other_face;
        e_to_e[other_cell][other_face] = cell_index;
        e_to_f[other_cell][other_face] = face_index;
      }
    }
  }

  m_e_to_e = std::move(e_to_e);
  m_e_to_f = std::move(e_to_f);
}
