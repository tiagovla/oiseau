// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <stdexcept>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/mesh/cell.hpp"

namespace oiseau::dg {

/**
 * @brief Returns the nodal reference element type for a mesh cell kind.
 * @param kind Mesh cell kind to convert.
 * @return Matching nodal reference element type.
 * @throws std::invalid_argument If the cell kind is unsupported.
 */
inline nodal::RefElementType ref_type_for(mesh::CellKind kind) {
  switch (kind) {
  case mesh::CellKind::Interval:
    return nodal::RefElementType::Line;
  case mesh::CellKind::Triangle:
    return nodal::RefElementType::Triangle;
  case mesh::CellKind::Quadrilateral:
    return nodal::RefElementType::Quadrilateral;
  case mesh::CellKind::Tetrahedron:
    return nodal::RefElementType::Tetrahedron;
  case mesh::CellKind::Hexahedron:
    return nodal::RefElementType::Hexahedron;
  default:
    throw std::invalid_argument("Unsupported mesh cell kind for DG reference element");
  }
}

}  // namespace oiseau::dg
