// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "oiseau/dg/dg_integration.hpp"
#include "oiseau/mesh/cell.hpp"

namespace oiseau::dg::detail {

DGIntegration::QuadratureRule make_face_quadrature(mesh::CellKind face_kind, unsigned order);

}  // namespace oiseau::dg::detail
