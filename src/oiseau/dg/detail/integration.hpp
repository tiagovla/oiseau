// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <span>

#include "oiseau/dg/dg_integration.hpp"
#include "oiseau/mesh/cell.hpp"

namespace oiseau::dg::detail {

using FacePermutation = std::array<std::size_t, 4>;

DGIntegration::QuadratureRule make_face_quadrature(mesh::CellKind face_kind, unsigned order);
std::array<double, 4> corner_weights(std::size_t corners, const DGIntegration::Matrix& nodes,
                                     std::size_t q);
DGIntegration::Matrix face_points(std::span<const std::array<double, 3>> cell_vertices,
                                  std::span<const int> face_vertices, std::size_t dimension,
                                  const FacePermutation& permutation,
                                  const DGIntegration::QuadratureRule& quadrature);

}  // namespace oiseau::dg::detail
