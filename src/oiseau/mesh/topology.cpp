// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
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

namespace {

using FaceLocation = std::pair<std::size_t, std::size_t>;

std::vector<std::size_t> canonical_face_key(std::span<const int> local_vertices,
                                            std::span<const std::size_t> cell_vertices) {
  std::vector<std::size_t> key;
  key.reserve(local_vertices.size());
  for (const int local_vertex : local_vertices) {
    if (local_vertex < 0 || static_cast<std::size_t>(local_vertex) >= cell_vertices.size()) {
      throw std::invalid_argument("Cell connectivity does not match its cell type");
    }
    key.push_back(cell_vertices[static_cast<std::size_t>(local_vertex)]);
  }
  std::ranges::sort(key);
  return key;
}

void connect_faces(FaceLocation face, FaceLocation other,
                   std::vector<std::vector<std::size_t>>& element_to_element,
                   std::vector<std::vector<std::size_t>>& element_to_face) {
  const auto [cell_index, face_index] = face;
  const auto [other_cell, other_face] = other;
  if (other_cell == cell_index || element_to_element[other_cell][other_face] != other_cell) {
    throw std::invalid_argument("Non-manifold face is shared by more than two cells");
  }
  element_to_element[cell_index][face_index] = other_cell;
  element_to_face[cell_index][face_index] = other_face;
  element_to_element[other_cell][other_face] = cell_index;
  element_to_face[other_cell][other_face] = face_index;
}

}  // namespace

Topology::Topology() = default;
Topology::~Topology() = default;

Topology::Topology(std::vector<std::vector<std::size_t>>&& conn, std::vector<CellType>&& cell_types)
    : m_conn(std::move(conn)), m_cell_types(std::move(cell_types)) {};

std::span<CellType> Topology::cell_types() {
  m_dimension = -1;
  return m_cell_types;
};
std::span<const CellType> Topology::cell_types() const { return m_cell_types; };

std::span<std::vector<std::size_t>> Topology::conn() { return m_conn; };
std::span<const std::vector<std::size_t>> Topology::conn() const { return m_conn; };
std::span<std::vector<std::size_t>> Topology::e_to_e() { return m_e_to_e; };
std::span<std::vector<std::size_t>> Topology::e_to_f() { return m_e_to_f; };

std::size_t Topology::n_cells() const { return m_conn.size(); }

int Topology::dimension() const {
  if (m_dimension >= 0) return m_dimension;

  int dimension = 0;
  for (const auto* cell : m_cell_types) {
    if (cell == nullptr) throw std::invalid_argument("Cell type cannot be null");
    dimension = std::max(dimension, cell->dimension());
  }
  m_dimension = dimension;
  return m_dimension;
}

void Topology::calculate_connectivity() {
  if (m_conn.size() != m_cell_types.size()) {
    throw std::invalid_argument("Cell connectivity and cell type counts must match");
  }

  const auto mesh_dimension = dimension();
  if (mesh_dimension > 3) {
    throw std::invalid_argument("Face connectivity supports meshes of dimension at most three");
  }

  // Each canonical face key records the first cell/local-face that uses it.
  // A second use creates the two-way adjacency; a third is non-manifold.
  std::map<std::vector<std::size_t>, FaceLocation> face_owner;
  std::vector<std::vector<std::size_t>> element_to_element(m_conn.size());
  std::vector<std::vector<std::size_t>> element_to_face(m_conn.size());

  for (std::size_t cell_index = 0; cell_index < m_conn.size(); ++cell_index) {
    const auto* cell = m_cell_types[cell_index];
    // Gmsh meshes commonly include lower-dimensional boundary entities. They
    // are not volume cells and therefore do not participate in this adjacency.
    if (mesh_dimension == 0 || cell->dimension() != mesh_dimension) continue;

    const auto face_vertices = cell->get_entity_vertices(mesh_dimension - 1);
    const auto& cell_vertices = m_conn[cell_index];
    element_to_element[cell_index].resize(face_vertices.size(), cell_index);
    element_to_face[cell_index].resize(face_vertices.size());

    for (std::size_t face_index = 0; face_index < face_vertices.size(); ++face_index) {
      auto key = canonical_face_key(face_vertices[face_index], cell_vertices);
      element_to_face[cell_index][face_index] = face_index;

      const FaceLocation face{cell_index, face_index};
      const auto [owner, inserted] = face_owner.emplace(std::move(key), face);
      if (inserted) continue;
      connect_faces(face, owner->second, element_to_element, element_to_face);
    }
  }

  m_e_to_e = std::move(element_to_element);
  m_e_to_f = std::move(element_to_face);
}
