// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <stdexcept>

#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace {

using namespace oiseau::dg;
using namespace oiseau::dg::solver;
using namespace oiseau::mesh;

Mesh interval_mesh() {
  return {Topology({{0, 1}}, {get_cell_type(CellKind::Interval)}), Geometry({0, 2}, 1)};
}

TEST(test_linear_dg, validates_spatial_directions_and_flux_policy) {
  auto mesh = interval_mesh();
  const Field u{"u"};
  const EquationSystem unsupported({u}, {dt(u) == -dy(u)});
  const EquationSystem zero_term({u}, {dt(u) == -dx(u) + 0.0 * dy(u)});
  const auto speed = [](const auto&) { return 1.0; };

  EXPECT_THROW(LinearDG1D(mesh, {1}, unsupported, speed), std::invalid_argument);
  EXPECT_THROW(LinearDG1D(mesh, {1}, zero_term, LinearDG1D::NumericalFlux{}),
               std::invalid_argument);
  EXPECT_NO_THROW(LinearDG1D(mesh, {1}, zero_term, speed));
}

}  // namespace
