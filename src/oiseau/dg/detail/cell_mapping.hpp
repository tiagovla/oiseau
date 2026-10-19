// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "oiseau/mesh/cell.hpp"

namespace oiseau::dg::detail {

/// Reorders mesh connectivity to match the node order of the reference element.
inline std::vector<std::size_t> reference_ordered_vertices(
    mesh::CellType cell_type, std::span<const std::size_t> mesh_vertices) {
  if (cell_type == nullptr) throw std::invalid_argument("Mesh cell type cannot be null");
  if (mesh_vertices.size() != static_cast<std::size_t>(cell_type->num_sub_entities(0))) {
    throw std::invalid_argument("Cell connectivity does not match its cell type");
  }

  std::vector<std::size_t> vertices(mesh_vertices.begin(), mesh_vertices.end());
  const auto kind = cell_type->kind();
  if (kind == mesh::CellKind::Quadrilateral || kind == mesh::CellKind::Hexahedron) {
    // Tensor-product reference nodes vary in the first coordinate first;

    std::swap(vertices[2], vertices[3]);
    if (kind == mesh::CellKind::Hexahedron) std::swap(vertices[6], vertices[7]);
  }
  return vertices;
}

}  // namespace oiseau::dg::detail
