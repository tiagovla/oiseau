// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "oiseau/dg/dg_space.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace {

using oiseau::dg::DGSpace;
using oiseau::mesh::CellKind;
using oiseau::mesh::Geometry;
using oiseau::mesh::get_cell_type;
using oiseau::mesh::Mesh;
using oiseau::mesh::Topology;

class UnsupportedCell : public oiseau::mesh::Cell {
 public:
  UnsupportedCell() {
    m_kind = CellKind::Undefined;
    m_dim = 1;
  }
};

Mesh line_mesh() {
  // 0----1----2 (two interval cells; coordinates 2, 6, 12).
  return {Topology({{0, 1}, {1, 2}},
                   {get_cell_type(CellKind::Interval), get_cell_type(CellKind::Interval)}),
          Geometry({2.0, 6.0, 12.0}, 1)};
}

Mesh mixed_mesh() {
  // Quad, boundary edge (cell 1), then triangles {1,4,2} and {4,5,2}.
  // The boundary edge shares the quad bottom.
  // 3-------2-------5
  // |       |     / |
  // |       |   /   |
  // |       | /     |
  // 0-------1-------4
  return {Topology({{0, 1, 2, 3}, {0, 1}, {1, 4, 2}, {4, 5, 2}},
                   {get_cell_type(CellKind::Quadrilateral), get_cell_type(CellKind::Interval),
                    get_cell_type(CellKind::Triangle), get_cell_type(CellKind::Triangle)}),
          Geometry({0, 0, 1, 0, 1, 1, 0, 1, 2, 0, 2, 1}, 2)};
}

Mesh volume_mesh() {
  // Hex mesh vertices:
  //       6-------7
  //      /|      /|
  //     4-------5 |
  //     | 2-----|-3
  //     |/      |/
  //     0-------1
  // Tetrahedron attaches to hex face triangle {1,3,5}:
  //       5
  //      /|\
  //     3-+-8
  //      \|/
  //       1
  return {
      Topology({{0, 1, 3, 2, 4, 5, 7, 6}, {1, 8, 3, 5}},
               {get_cell_type(CellKind::Hexahedron), get_cell_type(CellKind::Tetrahedron)}),
      Geometry({0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 2, 0, 0},
               3)};
}

}  // namespace

TEST(test_dg_space, uniform_order_broadcasts_and_maps_order_two_line_nodes) {
  const auto mesh = line_mesh();
  const DGSpace space(mesh, 2);

  ASSERT_EQ(space.orders().size(), 2);
  EXPECT_EQ(space.orders()[0], 2);
  EXPECT_EQ(space.orders()[1], 2);
  const auto* group = space.element_group(CellKind::Interval, 2);
  ASSERT_NE(group, nullptr);
  ASSERT_EQ(group->nodes.dimension(), 3);
  ASSERT_EQ(group->nodes.shape()[0], 2);
  ASSERT_EQ(group->nodes.shape()[1], 3);
  ASSERT_EQ(group->nodes.shape()[2], 1);
  // RefLine order 2 is left endpoint, midpoint, right endpoint.
  std::array<std::array<double, 3>, 2> expected = {{{2.0, 4.0, 6.0}, {6.0, 9.0, 12.0}}};
  for (std::size_t cell = 0; cell < 2; ++cell) {
    for (std::size_t node = 0; node < 3; ++node) {
      EXPECT_NEAR(group->nodes(cell, node, 0), expected[cell][node], 1e-12);
    }
  }
}

TEST(test_dg_space, groups_by_kind_and_order_and_skips_interleaved_boundary) {
  const auto mesh = mixed_mesh();
  const DGSpace space(mesh, {2, 3, 2});

  ASSERT_EQ(space.orders().size(), 3);
  EXPECT_EQ(space.orders()[0], 2);
  EXPECT_EQ(space.orders()[1], 3);
  EXPECT_EQ(space.orders()[2], 2);
  ASSERT_EQ(space.element_groups().size(), 3);
  const auto* quad = space.element_group(CellKind::Quadrilateral, 2);
  ASSERT_NE(quad, nullptr);
  EXPECT_EQ(quad->cell_indices, (std::vector<std::size_t>{0}));
  ASSERT_EQ(quad->nodes.shape(), (xt::xarray<double>::shape_type{1, 9, 2}));

  const auto* order_three = space.element_group(CellKind::Triangle, 3);
  ASSERT_NE(order_three, nullptr);
  EXPECT_EQ(order_three->cell_indices, (std::vector<std::size_t>{2}));
  ASSERT_EQ(order_three->nodes.shape(), (xt::xarray<double>::shape_type{1, 10, 2}));

  const auto* order_two = space.element_group(CellKind::Triangle, 2);
  ASSERT_NE(order_two, nullptr);
  EXPECT_EQ(order_two->cell_indices, (std::vector<std::size_t>{3}));
  ASSERT_EQ(order_two->nodes.shape(), (xt::xarray<double>::shape_type{1, 6, 2}));
  EXPECT_EQ(space.element_group(CellKind::Triangle, 4), nullptr);
}

TEST(test_dg_space, maps_hex_and_tet_vertices_in_three_dimensions) {
  const auto mesh = volume_mesh();
  const DGSpace space(mesh, 1);

  ASSERT_EQ(space.element_groups().size(), 2);
  const auto expect_vertices = [&space](CellKind kind, const auto& expected) {
    const auto* group = space.element_group(kind, 1);
    ASSERT_NE(group, nullptr);
    ASSERT_EQ(group->nodes.shape(), (xt::xarray<double>::shape_type{1, expected.size(), 3}));
    for (std::size_t node = 0; node < expected.size(); ++node) {
      for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        EXPECT_NEAR(group->nodes(0, node, coordinate), expected[node][coordinate], 1e-12);
      }
    }
  };
  const std::array<std::array<double, 3>, 8> expected_hex = {
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}};
  const std::array<std::array<double, 3>, 4> expected_tet = {
      {{1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {1, 0, 1}}};
  expect_vertices(CellKind::Hexahedron, expected_hex);
  expect_vertices(CellKind::Tetrahedron, expected_tet);
}

TEST(test_dg_space, rejects_invalid_orders_cell_types_and_vertices) {
  const auto mesh = mixed_mesh();
  const std::vector<unsigned> missing = {2, 3};
  const std::vector<unsigned> includes_boundary = {2, 3, 2, 1};
  const std::vector<unsigned> zero_on_later_cell = {2, 0, 2};
  EXPECT_THROW(DGSpace(mesh, missing), std::invalid_argument);
  EXPECT_THROW(DGSpace(mesh, includes_boundary), std::invalid_argument);
  EXPECT_THROW(DGSpace(mesh, zero_on_later_cell), std::invalid_argument);

  const UnsupportedCell unsupported_cell;
  const Mesh unsupported_mesh{Topology({{0, 1}}, {&unsupported_cell}), Geometry({0.0, 1.0}, 1)};
  EXPECT_THROW(DGSpace(unsupported_mesh, 1), std::invalid_argument);

  // Connectivity references vertex 1, but geometry contains only vertex 0.
  const Mesh invalid_vertex{Topology({{0, 1}}, {get_cell_type(CellKind::Interval)}),
                            Geometry({0.0}, 1)};
  EXPECT_THROW(DGSpace(invalid_vertex, 1), std::out_of_range);
}
