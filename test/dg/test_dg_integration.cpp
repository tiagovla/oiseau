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
#include <numeric>
#include <vector>

#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/dg_integration.hpp"
#include "oiseau/dg/dg_space.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace {

using namespace oiseau::dg;

// One quadrilateral and two triangles, followed by a boundary edge.
// Apply a shear and translation so the physical Jacobian is not diagonal.
//
// 3-------2-------5
// |       |     / |
// |       |   /   |
// |       | /     |
// 0-------1-------4

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
  for (const auto& point :
       std::array<std::array<double, 2>, 6>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {2, 1}}}) {
    coordinates.push_back(1.0 + 2.0 * point[0] + 0.5 * point[1]);
    coordinates.push_back(-2.0 + 0.25 * point[0] + 1.5 * point[1]);
  }
  return {Topology(std::move(cells), {quad, triangle, triangle, get_cell_type(CellKind::Interval)}),
          Geometry(std::move(coordinates), 2)};
}

class DGIntegrationTest : public ::testing::TestWithParam<bool> {};

TEST_P(DGIntegrationTest, mass_and_stiffness_satisfy_physical_integral_identities) {
  auto mesh = mixed_mesh(GetParam());
  const DGSpace space(mesh, {2, 3, 4});
  const DGConnectivity connectivity(mesh);
  const DGGeometry geometry(mesh);
  const DGIntegration integration(mesh, space, connectivity, geometry);

  ASSERT_EQ(integration.cells().size(), 4);
  EXPECT_EQ(integration.cell(3).group_index, DGConnectivity::no_neighbor);
  EXPECT_TRUE(integration.cell(3).faces.empty());
  for (std::size_t cell_index = 0; cell_index < 3; ++cell_index) {
    SCOPED_TRACE(cell_index);
    const auto& cell = integration.cell(cell_index);
    const auto mass = cell.mass();
    const auto inverse_mass = cell.inverse_mass();
    const auto stiffness_x = cell.stiffness_x();
    const auto stiffness_y = cell.stiffness_y();
    const auto nodes = cell.node_count();
    const auto area = cell_index == 0 ? 2.875 : 1.4375;
    EXPECT_NEAR(std::accumulate(mass.begin(), mass.end(), 0.0), area, 1e-11);
    for (std::size_t i = 0; i < nodes; ++i) {
      for (std::size_t j = 0; j < nodes; ++j) {
        double identity = 0.0;
        std::array<double, 2> boundary{};
        for (std::size_t k = 0; k < nodes; ++k) {
          identity += mass(i, k) * inverse_mass(k, j);
        }
        for (const auto& face : cell.faces) {
          for (std::size_t q = 0; q < face.point_count(); ++q) {
            const auto value =
                face.minus->trace(q, i) * face.physical_weight(q) * face.minus->trace(q, j);
            for (std::size_t d = 0; d < 2; ++d) boundary[d] += face.outward_normal[d] * value;
          }
        }
        EXPECT_NEAR(identity, i == j ? 1.0 : 0.0, 1e-11);
        EXPECT_NEAR(stiffness_x(i, j) + stiffness_x(j, i), boundary[0], 1e-11);
        EXPECT_NEAR(stiffness_y(i, j) + stiffness_y(j, i), boundary[1], 1e-11);
      }
    }
  }
}

TEST_P(DGIntegrationTest, unequal_orders_trace_the_same_polynomial_at_shared_face_points) {
  auto mesh = mixed_mesh(GetParam());
  const DGSpace space(mesh, {2, 3, 4});
  const DGConnectivity connectivity(mesh);
  const DGGeometry geometry(mesh);
  const DGIntegration integration(mesh, space, connectivity, geometry);
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
  auto mesh = mixed_mesh(GetParam());
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
  const auto integration = [this] {
    auto mesh = mixed_mesh(GetParam());
    const DGSpace space(mesh, 3);
    const DGConnectivity connectivity(mesh);
    const DGGeometry geometry(mesh);
    return DGIntegration(mesh, space, connectivity, geometry);
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
  auto mesh = mixed_mesh(GetParam());
  const DGSpace space(mesh, 3);
  const DGConnectivity connectivity(mesh);
  const DGGeometry geometry(mesh);
  const DGIntegration integration(mesh, space, connectivity, geometry);
  const auto quadrature = integration.cell(0).faces[0].minus->quadrature;
  for (std::size_t face = 0; face < 3; ++face) {
    EXPECT_EQ(integration.cell(1).faces[face].minus, integration.cell(2).faces[face].minus);
  }
  bool forward = false, reversed = false;
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

INSTANTIATE_TEST_SUITE_P(Orientations, DGIntegrationTest, ::testing::Bool());

}  // namespace
