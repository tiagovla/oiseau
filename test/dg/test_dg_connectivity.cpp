// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace {

using oiseau::dg::DGConnectivity;
using oiseau::mesh::CellKind;
using oiseau::mesh::Geometry;
using oiseau::mesh::get_cell_type;
using oiseau::mesh::Mesh;
using oiseau::mesh::Topology;

Mesh mixed_mesh(std::vector<std::size_t> triangle_vertices) {
  // Cell 0 is the quad. Cell 1 is its bottom edge. Cell 2 is the triangle.
  // 3-------2
  // |       | \
  // |       |  4
  // |       | /
  // 0-------1
  return {Topology({{0, 1, 2, 3}, {0, 1}, std::move(triangle_vertices)},
                   {get_cell_type(CellKind::Quadrilateral), get_cell_type(CellKind::Interval),
                    get_cell_type(CellKind::Triangle)}),
          Geometry({0, 0, 1, 0, 1, 2, 0, 2, 2, 1}, 2)};
}

}  // namespace

TEST(test_dg_connectivity, mixed_cells_store_both_directions_of_a_reversed_shared_edge) {
  auto mesh = mixed_mesh({2, 4, 1});
  const DGConnectivity connectivity(mesh);

  const auto quad_faces = connectivity.faces(0);
  const auto triangle_faces = connectivity.faces(2);
  ASSERT_EQ(quad_faces.size(), 4);
  ASSERT_EQ(triangle_faces.size(), 3);

  // Local face 1 runs from vertex 1 to 2 on the quad and from 2 to 1 on the triangle.
  const auto& quad_to_triangle = quad_faces[1];
  EXPECT_FALSE(quad_to_triangle.is_boundary());
  EXPECT_EQ(quad_to_triangle.neighbor_cell, 2);
  EXPECT_EQ(quad_to_triangle.neighbor_face, 1);
  EXPECT_EQ(quad_to_triangle.vertex_permutation, (std::vector<std::size_t>{1, 0}));

  const auto& triangle_to_quad = triangle_faces[1];
  EXPECT_FALSE(triangle_to_quad.is_boundary());
  EXPECT_EQ(triangle_to_quad.neighbor_cell, 0);
  EXPECT_EQ(triangle_to_quad.neighbor_face, 1);
  EXPECT_EQ(triangle_to_quad.vertex_permutation, (std::vector<std::size_t>{1, 0}));
}

TEST(test_dg_connectivity, aligned_shared_edges_keep_the_vertex_order) {
  auto mesh = mixed_mesh({4, 1, 2});
  const DGConnectivity connectivity(mesh);

  // The quad's face 1 and the triangle's face 0 both list vertices {1, 2}.
  const auto quad_faces = connectivity.faces(0);
  const auto triangle_faces = connectivity.faces(2);
  ASSERT_EQ(quad_faces.size(), 4);
  ASSERT_EQ(triangle_faces.size(), 3);
  EXPECT_EQ(quad_faces[1].neighbor_cell, 2);
  EXPECT_EQ(quad_faces[1].neighbor_face, 0);
  EXPECT_EQ(triangle_faces[0].neighbor_cell, 0);
  EXPECT_EQ(triangle_faces[0].neighbor_face, 1);
  const std::vector<std::size_t> same_order{0, 1};
  EXPECT_EQ(quad_faces[1].vertex_permutation, same_order);
  EXPECT_EQ(triangle_faces[0].vertex_permutation, same_order);
}

TEST(test_dg_connectivity, a_boundary_face_has_a_record_and_the_mesh_edge_has_an_empty_span) {
  auto mesh = mixed_mesh({2, 4, 1});
  const DGConnectivity connectivity(mesh);

  const auto quad_faces = connectivity.faces(0);
  ASSERT_EQ(quad_faces.size(), 4);
  const auto& bottom_face = quad_faces[0];
  EXPECT_TRUE(bottom_face.is_boundary());
  EXPECT_EQ(bottom_face.neighbor_cell, DGConnectivity::no_neighbor);
  EXPECT_EQ(bottom_face.neighbor_face, DGConnectivity::no_neighbor);
  EXPECT_TRUE(bottom_face.vertex_permutation.empty());
  EXPECT_TRUE(connectivity.faces(1).empty());
}

TEST(test_dg_connectivity, face_lookup_rejects_an_out_of_range_cell_index) {
  auto mesh = mixed_mesh({2, 4, 1});
  const DGConnectivity connectivity(mesh);

  EXPECT_THROW(connectivity.faces(3), std::out_of_range);
}

TEST(test_dg_connectivity, rejects_a_one_dimensional_mesh) {
  Mesh mesh{Topology({{0, 1}}, {get_cell_type(CellKind::Interval)}), Geometry({0, 1}, 1)};

  EXPECT_THROW(DGConnectivity{mesh}, std::invalid_argument);
}

TEST(test_dg_connectivity, rejects_an_edge_shared_by_three_cells) {
  const auto triangle = get_cell_type(CellKind::Triangle);
  // Three triangles all use the edge between vertices 0 and 1.
  Mesh mesh{Topology({{0, 1, 2}, {1, 0, 3}, {0, 1, 4}}, {triangle, triangle, triangle}),
            Geometry({0, 0, 2, 0, 0, 1, 0, -1, 2, 1}, 2)};

  EXPECT_THROW(DGConnectivity{mesh}, std::invalid_argument);
}
