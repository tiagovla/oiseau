// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/dg/solver/state.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"
#include "oiseau/utils/integration.hpp"

namespace oiseau {
namespace {

using namespace oiseau::dg;
using namespace oiseau::dg::solver;
using namespace oiseau::mesh;

Mesh interval_mesh() {
  return {Topology({{0, 1}}, {get_cell_type(CellKind::Interval)}),
          oiseau::mesh::Geometry({0, 2}, 1)};
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

  EXPECT_THROW(solver.rhs(State{}), std::invalid_argument);
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
  const double* boundary_buffer = nullptr;
  const LinearDG1D solver(mesh, {1}, equations, flux,
                          [&boundary_buffer](auto interior, const auto&, auto exterior) {
                            if (!boundary_buffer) boundary_buffer = exterior.data();
                            EXPECT_EQ(exterior.data(), boundary_buffer);
                            EXPECT_DOUBLE_EQ(exterior[0], interior[0]);
                            EXPECT_DOUBLE_EQ(exterior[1], interior[1]);
                            exterior[0] += 1.0;
                            exterior[1] -= 2.0;
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

TEST(test_linear_dg, rejects_flux_with_the_wrong_number_of_fields) {
  auto mesh = interval_mesh();
  const Field u{"u"};
  const EquationSystem equations({u}, {dt(u) == -dx(u)});
  const auto empty_result = [](const auto&...) { return std::vector<double>{}; };
  const LinearDG1D bad_flux(mesh, {1}, equations, LinearDG1D::NumericalFlux(empty_result));

  EXPECT_THROW(bad_flux.rhs(bad_flux.zero_state()), std::invalid_argument);
}

TEST(test_linear_dg, differentiates_a_linear_field_across_rotated_3d_cells) {
  const auto tetrahedron = get_cell_type(CellKind::Tetrahedron);
  Mesh mesh{Topology({{0, 1, 2, 3}, {2, 0, 1, 4}}, {tetrahedron, tetrahedron}),
            oiseau::mesh::Geometry({0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1}, 3)};
  const Field u{"u"};
  const EquationSystem equations({u}, {dt(u) == -dx(u) - 2.0 * dy(u) - 3.0 * dz(u)});
  const LinearDG3D solver(mesh, {1, 2}, equations, [](const auto&) { return 4.0; });
  const auto state =
      solver.initial_state(u, [](double x, double y, double z) { return x + 2.0 * y + 3.0 * z; });

  const auto rhs = solver.rhs(state);
  ASSERT_EQ(rhs.groups.size(), 2);
  for (const auto& group : rhs.groups) {
    for (const auto value : group) EXPECT_NEAR(value, -14.0, 1e-11);
  }
}

TEST(test_linear_dg, quadratic_advection_converges_at_third_order_in_space) {
  constexpr unsigned order = 2;
  constexpr double final_time = 0.125;
  constexpr std::size_t steps = 256;
  const Field u{"u"};
  const EquationSystem equations({u}, {dt(u) == -dx(u)});
  const auto pulse = [](double x) {
    // The smooth pulse stays inside [0, 1], so the exact inflow is zero.
    const double scaled = (x - 0.375) / 0.25;
    if (std::abs(scaled) >= 1.0) return 0.0;
    return std::exp(1.0 - 1.0 / (1.0 - scaled * scaled));
  };

  const auto [points, weights] = oiseau::utils::integration::quadrature(15);
  const auto reference = nodal::get_ref_element(nodal::RefElementType::Line, order);
  const auto interpolation = xt::linalg::dot(reference->vandermonde(points), reference->inv_v());
  const auto l2_error = [&](std::size_t cell_count) {
    const double h = 1.0 / static_cast<double>(cell_count);
    std::vector<double> vertices(cell_count + 1);
    std::vector<std::vector<std::size_t>> cells(cell_count);
    for (std::size_t i = 0; i < cell_count; ++i) {
      vertices[i] = static_cast<double>(i) * h;
      cells[i] = {i, i + 1};
    }
    vertices.back() = 1.0;
    const auto interval = get_cell_type(CellKind::Interval);
    Mesh mesh{Topology(std::move(cells), std::vector<CellType>(cell_count, interval)),
              oiseau::mesh::Geometry(std::move(vertices), 1)};
    const Space space(mesh, order);
    const LinearDG1D solver(
        mesh, space, equations, [](const auto&) { return 1.0; },
        [](const auto&, const auto&, auto exterior) { exterior[0] = 0.0; });
    auto state = solver.initial_state(u, pulse);
    oiseau::time::RK4 rk4;
    oiseau::time::run(rk4, solver, state, 0.0, final_time / steps, steps);

    const auto groups = space.element_groups();
    const auto& cell_indices = groups.front().cell_indices;
    const auto& values = state.groups.front();
    double squared_error = 0.0;
    for (std::size_t cell = 0; cell < cell_count; ++cell) {
      const double left = static_cast<double>(cell_indices[cell]) * h;
      const auto interpolated =
          xt::linalg::dot(interpolation, xt::view(values, cell, xt::all(), 0));
      for (std::size_t q = 0; q < weights.size(); ++q) {
        const double x = left + 0.5 * h * (1.0 + points(q));
        const double error = interpolated(q) - pulse(x - final_time);
        squared_error += 0.5 * h * weights(q) * error * error;
      }
    }
    return std::sqrt(squared_error);
  };

  // Use the same small RK4 step on every mesh to keep time error below spatial error.
  double coarse_error = l2_error(32);
  for (const std::size_t cell_count : {64, 128}) {
    const double fine_error = l2_error(cell_count);
    ASSERT_TRUE(std::isfinite(fine_error));
    ASSERT_GT(fine_error, 0.0);
    EXPECT_GE(std::log2(coarse_error / fine_error), 2.75)
        << "at " << cell_count << " cells, L2 error " << coarse_error << " -> " << fine_error;
    coarse_error = fine_error;
  }
}

}  // namespace
}  // namespace oiseau
