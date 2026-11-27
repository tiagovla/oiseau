// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xoperation.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/connectivity.hpp"
#include "oiseau/dg/detail/face_quadrature.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/geometry.hpp"
#include "oiseau/dg/integration.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau {
namespace {

using namespace oiseau::dg;
using oiseau::mesh::CellKind;

TEST(DGFaceQuadratureTest, point_has_one_weight_and_no_coordinates) {
  const auto [nodes, weights] = detail::make_face_quadrature(CellKind::Point, 2);
  ASSERT_EQ(nodes.dimension(), 2);
  EXPECT_EQ(nodes.shape()[0], 1);
  EXPECT_EQ(nodes.shape()[1], 0);
  ASSERT_EQ(weights.size(), 1);
  EXPECT_DOUBLE_EQ(weights(0), 1.0);
}

TEST(DGFaceQuadratureTest, line_integrates_fourth_degree_moment) {
  const auto [nodes, weights] = detail::make_face_quadrature(CellKind::Interval, 2);
  ASSERT_EQ(nodes.dimension(), 1);
  ASSERT_EQ(nodes.size(), weights.size());
  EXPECT_NEAR(xt::sum(weights)(), 2.0, 1e-12);
  EXPECT_NEAR(xt::sum(xt::pow(nodes, 4) * weights)(), 2.0 / 5.0, 1e-12);
}

TEST(DGFaceQuadratureTest, triangle_integrates_mixed_fifth_degree_moment) {
  const auto [nodes, weights] = detail::make_face_quadrature(CellKind::Triangle, 2);
  ASSERT_EQ(nodes.dimension(), 2);
  ASSERT_EQ(nodes.shape()[0], weights.size());
  ASSERT_EQ(nodes.shape()[1], 2);
  const auto r = xt::col(nodes, 0);
  const auto s = xt::col(nodes, 1);
  EXPECT_NEAR(xt::sum(weights)(), 2.0, 1e-12);
  // Integrate r^3 s^2 over -1 <= r <= 1 and -1 <= s <= -r.
  EXPECT_NEAR(xt::sum(xt::pow(r, 3) * xt::pow(s, 2) * weights)(), -2.0 / 21.0, 1e-12);
}

TEST(DGFaceQuadratureTest, square_integrates_tensor_product_moment) {
  const auto [nodes, weights] = detail::make_face_quadrature(CellKind::Quadrilateral, 2);
  ASSERT_EQ(nodes.dimension(), 2);
  ASSERT_EQ(nodes.shape()[0], weights.size());
  ASSERT_EQ(nodes.shape()[1], 2);
  const auto r = xt::col(nodes, 0);
  const auto s = xt::col(nodes, 1);
  EXPECT_NEAR(xt::sum(weights)(), 4.0, 1e-12);
  EXPECT_NEAR(xt::sum(xt::pow(r, 4) * xt::pow(s, 4) * weights)(), 4.0 / 25.0, 1e-12);
}

TEST(DGFaceQuadratureTest, rejects_invalid_shape_and_order) {
  EXPECT_THROW(detail::make_face_quadrature(CellKind::Tetrahedron, 2), std::logic_error);
  const auto oversized_order = std::numeric_limits<unsigned>::max();
  EXPECT_THROW(detail::make_face_quadrature(CellKind::Interval, oversized_order),
               std::invalid_argument);
}

TEST(DGFaceQuadratureTest, corner_weights_map_a_point_on_an_oriented_edge) {
  const std::array cell_vertices = {std::array{-1.0, -1.0, 0.0}, std::array{1.0, -1.0, 0.0},
                                    std::array{-1.0, 1.0, 0.0}};
  const std::array face_vertices{1, 2};
  const Integration::QuadratureRule quadrature{.nodes = {0.5}, .weights = {1.0}};

  EXPECT_EQ(detail::corner_weights(face_vertices.size(), quadrature.nodes, 0),
            (std::array{0.25, 0.75, 0.0, 0.0}));
  const auto points =
      detail::face_points(cell_vertices, face_vertices, 2, {0, 1, 2, 3}, quadrature);
  const auto reversed =
      detail::face_points(cell_vertices, face_vertices, 2, {1, 0, 2, 3}, quadrature);
  EXPECT_EQ(points, (Integration::Matrix{{-0.5, 0.5}}));
  EXPECT_EQ(reversed, (Integration::Matrix{{0.5, -0.5}}));
}

// One quadrilateral and two triangles, followed by a boundary edge.
// Apply a shear and translation so the physical Jacobian is not diagonal.
//
// 3-------2-------5
// |       | \     |
// |       |   \   |
// |       |     \ |
// 0-------1-------4
//
// Vertex orders are (0, 1, 2, 3), (1, 4, 2), and (4, 5, 2).
//
// Cell winding for reversed=true
//
// 3-------2-------5
// |       | \ CCW |
// |  CW   |   \   |
// |       | CW  \ |
// 0-------1-------4
//
// Vertex orders are (3, 2, 1, 0), (2, 4, 1), and (4, 5, 2).

oiseau::mesh::Mesh mixed_mesh(bool reversed) {
  using namespace oiseau::mesh;
  std::vector<std::vector<std::size_t>> cells = {{0, 1, 2, 3}, {1, 4, 2}, {4, 5, 2}, {0, 1}};
  if (reversed) {
    // Leave the second triangle positive to exercise both face permutations.
    std::ranges::reverse(cells[0]);
    std::ranges::reverse(cells[1]);
  }
  const auto quad = get_cell_type(CellKind::Quadrilateral);
  const auto triangle = get_cell_type(CellKind::Triangle);
  std::vector<double> coordinates;
  // Affine transformation: x' = 1 + 2x + 0.5y, y' = -2 + 0.25x + 1.5y.
  for (const auto& point :
       std::array<std::array<double, 2>, 6>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {2, 1}}}) {
    coordinates.push_back(1.0 + 2.0 * point[0] + 0.5 * point[1]);
    coordinates.push_back(-2.0 + 0.25 * point[0] + 1.5 * point[1]);
  }
  return {Topology(std::move(cells), {quad, triangle, triangle, get_cell_type(CellKind::Interval)}),
          oiseau::mesh::Geometry(std::move(coordinates), 2)};
}

class DGIntegrationTest : public ::testing::TestWithParam<bool> {};

TEST_P(DGIntegrationTest, mass_and_stiffness_satisfy_physical_integral_identities) {
  bool reversed = GetParam();
  auto mesh = mixed_mesh(reversed);
  const Space space(mesh, {2, 3, 4});
  const Connectivity connectivity(mesh);
  const Geometry geometry(mesh);
  const Integration integration(mesh, space, connectivity, geometry);

  constexpr double tolerance = 1e-11;
  const auto expect_matrix_near = [](const char* name, const Integration::Matrix& actual,
                                     const Integration::Matrix& expected) {
    SCOPED_TRACE(name);
    ASSERT_EQ(actual.shape(), expected.shape());
    for (const auto& index : xt::argwhere(!xt::isclose(actual, expected, 0.0, tolerance))) {
      const auto i = index[0];
      const auto j = index[1];
      EXPECT_NEAR(actual(i, j), expected(i, j), tolerance) << "entry (" << i << ", " << j << ")";
    }
  };

  ASSERT_EQ(integration.cells().size(), 4);
  for (const auto [cell_index, cell] : std::views::enumerate(integration.cells())) {
    SCOPED_TRACE(cell_index);
    ASSERT_EQ(cell.has_group(), cell_index != 3);
    if (!cell.has_group()) {
      EXPECT_TRUE(cell.faces.empty());
      continue;
    }
    const auto mass = cell.mass();
    const auto inverse_mass = cell.inverse_mass();
    const auto stiffness_x = cell.stiffness_x();
    const auto stiffness_y = cell.stiffness_y();
    const auto nodes = cell.node_count();
    const auto area = cell_index == 0 ? 2.875 : 1.4375;
    EXPECT_NEAR(xt::sum(mass)(), area, tolerance);
    expect_matrix_near("mass times inverse mass", xt::linalg::dot(mass, inverse_mass),
                       xt::eye<double>(nodes));

    Integration::Matrix boundary_x = xt::zeros<double>({nodes, nodes});
    Integration::Matrix boundary_y = xt::zeros<double>({nodes, nodes});
    for (const auto& face : cell.faces) {
      const auto& trace = face.minus->trace;
      auto weights = xt::empty<double>({face.point_count()});
      for (std::size_t q = 0; q < face.point_count(); ++q) {
        weights(q) = face.physical_weight(q);
      }
      const auto weighted_trace = xt::view(weights, xt::all(), xt::newaxis()) * trace;
      const auto face_mass = xt::linalg::dot(xt::transpose(trace), weighted_trace);
      boundary_x += face.outward_normal[0] * face_mass;
      boundary_y += face.outward_normal[1] * face_mass;
    }
    expect_matrix_near("x stiffness", stiffness_x + xt::transpose(stiffness_x), boundary_x);
    expect_matrix_near("y stiffness", stiffness_y + xt::transpose(stiffness_y), boundary_y);
  }
}

TEST_P(DGIntegrationTest, unequal_orders_trace_the_same_polynomial_at_shared_face_points) {
  bool reversed = GetParam();
  auto mesh = mixed_mesh(reversed);
  const Space space(mesh, {2, 3, 4});
  const Connectivity connectivity(mesh);
  const Geometry geometry(mesh);
  const Integration integration(mesh, space, connectivity, geometry);
  const auto groups = space.element_groups();
  std::size_t interior_sides = 0;
  for (std::size_t cell_index = 0; cell_index < 3; ++cell_index) {
    const auto& cell = integration.cell(cell_index);
    for (const auto& face : cell.faces) {
      if (face.is_boundary()) continue;
      ++interior_sides;
      const auto& neighbor = integration.cell(face.neighbor_cell);
      const auto& minus_nodes = groups[cell.group_index].nodes;
      const auto& plus_nodes = groups[neighbor.group_index].nodes;
      ASSERT_EQ(face.minus->trace.shape()[0], face.plus->trace.shape()[0]);
      for (std::size_t q = 0; q < face.point_count(); ++q) {
        // Constants, both coordinates, and a quadratic distinguish alignment errors.
        for (int polynomial = 0; polynomial < 4; ++polynomial) {
          const auto value = [polynomial](double x, double y) {
            if (polynomial == 0) return 1.0;
            if (polynomial == 1) return x;
            if (polynomial == 2) return y;
            return x * x + x * y + y * y;
          };
          double minus = 0.0, plus = 0.0;
          for (std::size_t i = 0; i < minus_nodes.shape()[1]; ++i) {
            minus += face.minus->trace(q, i) * value(minus_nodes(cell.group_cell_index, i, 0),
                                                     minus_nodes(cell.group_cell_index, i, 1));
          }
          for (std::size_t i = 0; i < plus_nodes.shape()[1]; ++i) {
            plus += face.plus->trace(q, i) * value(plus_nodes(neighbor.group_cell_index, i, 0),
                                                   plus_nodes(neighbor.group_cell_index, i, 1));
          }
          EXPECT_NEAR(minus, plus, 1e-11);
        }
      }
    }
  }
  EXPECT_EQ(interior_sides, 4);
}

TEST_P(DGIntegrationTest, linear_solver_preserves_constants_and_differentiates_quadratics) {
  bool reversed = GetParam();
  auto mesh = mixed_mesh(reversed);
  const Field u{"u"};
  const EquationSystem equations({u}, {dt(u) == -dx(u) - 0.5 * dy(u)});
  const solver::LinearDG2D solver(mesh, {2, 3, 4}, equations, [](const auto& normal) {
    return std::abs(normal[0] + 0.5 * normal[1]);
  });
  for (const bool constant : {true, false}) {
    SCOPED_TRACE(constant);
    const auto state = solver.initial_state(
        u, [constant](double x, double y) { return constant ? 3.0 : x * x + x * y + y * y; });
    const auto rhs = solver.rhs(state);
    const auto groups = solver.space().element_groups();
    for (std::size_t g = 0; g < groups.size(); ++g) {
      const auto& nodes = groups[g].nodes;
      for (std::size_t cell = 0; cell < nodes.shape()[0]; ++cell) {
        for (std::size_t node = 0; node < nodes.shape()[1]; ++node) {
          const auto x = nodes(cell, node, 0), y = nodes(cell, node, 1);
          EXPECT_NEAR(rhs.groups[g](cell, node, 0), constant ? 0.0 : -2.5 * x - 2.0 * y, 1e-10);
        }
      }
    }
  }
}

TEST_P(DGIntegrationTest, reference_operators_are_shared_and_own_their_lifetime) {
  const auto integration = [] {
    bool reversed = GetParam();
    auto mesh = mixed_mesh(reversed);
    const Space space(mesh, 3);
    return Integration(mesh, space);
  }();
  const auto copy = integration;
  EXPECT_EQ(integration.cell(1).reference, integration.cell(2).reference);
  EXPECT_NE(integration.cell(0).reference, integration.cell(1).reference);
  EXPECT_EQ(copy.cell(1).reference, integration.cell(1).reference);
  EXPECT_EQ(copy.cell(1).faces[0].minus, integration.cell(1).faces[0].minus);
  EXPECT_EQ(copy.cell(1).faces[0].point_count(), 4);
  const auto mass = copy.cell(1).mass();
  EXPECT_NEAR(std::accumulate(mass.begin(), mass.end(), 0.0), 1.4375, 1e-11);
}

TEST_P(DGIntegrationTest, face_rules_share_reference_data_and_apply_cell_scaling) {
  bool reversed = GetParam();
  auto mesh = mixed_mesh(reversed);
  const Space space(mesh, 3);
  const Connectivity connectivity(mesh);
  const Geometry geometry(mesh);
  const Integration integration(mesh, space, connectivity, geometry);
  const auto quadrature = integration.cell(0).faces[0].minus->quadrature;
  for (std::size_t face = 0; face < 3; ++face) {
    EXPECT_EQ(integration.cell(1).faces[face].minus, integration.cell(2).faces[face].minus);
  }
  bool forward = false;
  reversed = false;
  for (std::size_t cell_index = 0; cell_index < 3; ++cell_index) {
    const auto& cell = integration.cell(cell_index);
    const auto inverse_mass = cell.inverse_mass();
    for (std::size_t f = 0; f < cell.faces.size(); ++f) {
      const auto& face = cell.faces[f];
      EXPECT_EQ(face.minus->quadrature, quadrature);
      EXPECT_EQ(face.point_count(), 4);
      double measure = 0.0;
      for (std::size_t q = 0; q < face.point_count(); ++q) {
        EXPECT_GT(face.physical_weight(q), 0.0);
        measure += face.physical_weight(q);
        for (std::size_t i = 0; i < cell.node_count(); ++i) {
          double lift = 0.0;
          for (std::size_t j = 0; j < cell.node_count(); ++j) {
            lift += inverse_mass(i, j) * face.minus->trace(q, j) * face.physical_weight(q);
          }
          EXPECT_NEAR(face.physical_lift(q, i), lift, 1e-11);
        }
      }
      EXPECT_NEAR(measure, geometry.cell(cell_index).faces[f].measure, 1e-12);
      if (face.is_boundary()) {
        EXPECT_EQ(face.plus, nullptr);
        continue;
      }
      EXPECT_EQ(face.plus->quadrature, quadrature);
      const auto& connection = connectivity.faces(cell_index)[f];
      const auto& opposite = integration.cell(face.neighbor_cell).faces[connection.neighbor_face];
      if (connection.vertex_permutation[0] == 0) {
        forward = true;
        EXPECT_EQ(face.plus, opposite.minus);
      } else {
        reversed = true;
        EXPECT_NE(face.plus, opposite.minus);
        for (std::size_t q = 0; q < face.point_count(); ++q) {
          for (std::size_t node = 0; node < face.plus->trace.shape()[1]; ++node) {
            EXPECT_NEAR(face.plus->trace(q, node),
                        opposite.minus->trace(face.point_count() - 1 - q, node), 1e-12);
          }
        }
      }
    }
  }
  EXPECT_TRUE(forward);
  EXPECT_EQ(reversed, GetParam());
}

TEST(DGIntegrationInputTest, rejects_inputs_from_different_meshes) {
  auto mesh = mixed_mesh(false);
  auto other_mesh = mesh;
  const Space space(mesh, 1);
  const Connectivity connectivity(mesh);
  const Geometry geometry(mesh);
  const Space other_space(other_mesh, 1);
  const Connectivity other_connectivity(other_mesh);
  const Geometry other_geometry(other_mesh);

  EXPECT_THROW(Integration(mesh, other_space, connectivity, geometry), std::invalid_argument);
  EXPECT_THROW(Integration(mesh, space, other_connectivity, geometry), std::invalid_argument);
  EXPECT_THROW(Integration(mesh, space, connectivity, other_geometry), std::invalid_argument);
}

TEST(DGIntegrationDimensionTest, interval_has_only_x_volume_operators) {
  using namespace oiseau::mesh;
  const auto interval = get_cell_type(CellKind::Interval);
  Mesh mesh{Topology({{0, 1}}, {interval}), oiseau::mesh::Geometry({0.0, 2.0}, 1)};
  const Space space(mesh, 1);
  const Connectivity connectivity(mesh);
  const oiseau::dg::Geometry geometry(mesh);
  const Integration integration(mesh, space, connectivity, geometry);
  const auto& cell = integration.cell(0);

  EXPECT_EQ(cell.weak_derivative_x.size(), 4);
  EXPECT_TRUE(cell.weak_derivative_y.empty());
  EXPECT_TRUE(cell.weak_derivative_z.empty());
  EXPECT_EQ(cell.stiffness_y().dimension(), 0);
  EXPECT_EQ(cell.stiffness_z().dimension(), 0);
  ASSERT_EQ(cell.faces.size(), 2);
  EXPECT_EQ(cell.faces[0].minus->quadrature, cell.faces[1].minus->quadrature);
}

TEST(DGIntegrationDimensionTest, rotated_tetrahedral_interface_shares_both_oriented_rules) {
  using namespace oiseau::mesh;
  const auto tetrahedron = get_cell_type(CellKind::Tetrahedron);
  Mesh mesh{Topology({{0, 1, 2, 3}, {2, 0, 1, 4}}, {tetrahedron, tetrahedron}),
            oiseau::mesh::Geometry({0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1}, 3)};
  const Space space(mesh, {1, 2});
  const Connectivity connectivity(mesh);
  const oiseau::dg::Geometry geometry(mesh);
  const Integration integration(mesh, space, connectivity, geometry);

  constexpr std::size_t shared_face = 3;
  const auto& connection = connectivity.faces(0)[shared_face];
  ASSERT_EQ(connection.neighbor_cell, 1);
  ASSERT_EQ(connection.neighbor_face, shared_face);
  EXPECT_EQ(connection.vertex_permutation, (std::vector<std::size_t>{1, 2, 0}));

  const auto& cell = integration.cell(0);
  const auto& first = cell.faces[shared_face];
  const auto& second = integration.cell(1).faces[shared_face];
  ASSERT_TRUE(first.plus);
  ASSERT_TRUE(second.plus);
  EXPECT_EQ(first.minus, second.plus);
  EXPECT_EQ(first.plus, second.minus);
  EXPECT_EQ(first.minus->quadrature, first.plus->quadrature);

  const auto weak_z = xt::linalg::dot(cell.inverse_mass(), cell.stiffness_z());
  ASSERT_EQ(cell.weak_derivative_z.size(), weak_z.size());
  for (const auto [i, expected] : std::views::enumerate(weak_z)) {
    EXPECT_NEAR(cell.weak_derivative_z[i], expected, 1e-12);
  }
}

INSTANTIATE_TEST_SUITE_P(Orientations, DGIntegrationTest, ::testing::Bool());

}  // namespace
}  // namespace oiseau
