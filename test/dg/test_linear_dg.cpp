// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>
#include <xtensor/containers/xarray.hpp>

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

TEST(test_linear_dg, validates_state_shape_and_strides_before_sampling_or_computing_rhs) {
  auto mesh = interval_mesh();
  const Field u{"u"}, v{"v"};
  const EquationSystem equations({u, v}, {dt(u) == -dx(v), dt(v) == -dx(u)});
  const LinearDG1D solver(mesh, {1}, equations, [](const auto&) { return 1.0; });

  EXPECT_THROW(solver.rhs(DGState{}), std::invalid_argument);
  auto state = solver.zero_state();
  state.groups[0].reshape({2, 2});
  EXPECT_THROW(solver.rhs(state), std::invalid_argument);

  using Array = xt::xarray<double>;
  state.groups[0] = Array(Array::storage_type(4, 0.0), {1, 2, 2}, {0, 1, 2});
  EXPECT_THROW(solver.set_initial_condition(state, u, [](double x) { return x; }),
               std::invalid_argument);

  state = solver.initial_state(v, [](double x) { return x + 1.0; });
  EXPECT_THROW(solver.set_initial_condition(state, u, {}), std::invalid_argument);
  solver.set_initial_condition(state, u, [](double x) { return 2.0 * x; });
  EXPECT_DOUBLE_EQ(state.groups[0](0, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(state.groups[0](0, 1, 0), 4.0);
  EXPECT_DOUBLE_EQ(state.groups[0](0, 0, 1), 1.0);
  EXPECT_DOUBLE_EQ(state.groups[0](0, 1, 1), 3.0);
}

TEST(test_linear_dg, passes_both_normal_fluxes_to_the_policy_and_lifts_its_correction) {
  auto mesh = interval_mesh();
  const Field u{"u"}, v{"v"};
  const EquationSystem equations({u, v}, {dt(u) == -dx(v), dt(v) == -2.0 * dx(u)});
  std::vector<double> normals;
  const LinearDG1D::NumericalFlux flux =
      [&normals](std::span<const double> left, std::span<const double> right,
                 const LinearDG1D::Normal& normal, std::span<const double> left_flux,
                 std::span<const double> right_flux) {
        normals.push_back(normal[0]);
        EXPECT_NEAR(left[0], 2.0, 1e-12);
        EXPECT_NEAR(left[1], 3.0, 1e-12);
        EXPECT_NEAR(right[0], 3.0, 1e-12);
        EXPECT_NEAR(right[1], 1.0, 1e-12);
        EXPECT_NEAR(left_flux[0], 3.0 * normal[0], 1e-12);
        EXPECT_NEAR(left_flux[1], 4.0 * normal[0], 1e-12);
        EXPECT_NEAR(right_flux[0], normal[0], 1e-12);
        EXPECT_NEAR(right_flux[1], 6.0 * normal[0], 1e-12);
        return std::vector<double>{left_flux[0] + 1.0, left_flux[1] - 2.0};
      };
  const LinearDG1D solver(mesh, {1}, equations, flux, [](const auto& inside, const auto&) {
    return std::vector<double>{inside[0] + 1.0, inside[1] - 2.0};
  });
  auto state = solver.initial_state(u, [](double) { return 2.0; });
  solver.set_initial_condition(state, v, [](double) { return 3.0; });

  const auto rhs = solver.rhs(state);
  EXPECT_EQ(normals, (std::vector<double>{-1.0, 1.0}));
  for (std::size_t node = 0; node < 2; ++node) {
    EXPECT_NEAR(rhs.groups[0](0, node, 0), -1.0, 1e-12);
    EXPECT_NEAR(rhs.groups[0](0, node, 1), 2.0, 1e-12);
  }
}

}  // namespace
