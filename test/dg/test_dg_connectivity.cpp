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

#include "oiseau/dg/connectivity.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"

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

TEST(test_dg_connectivity, intervals_share_an_endpoint_and_ignore_point_tags) {
  const auto interval = get_cell_type(CellKind::Interval);
  const auto point = get_cell_type(CellKind::Point);
  Mesh mesh{Topology({{0, 1}, {0}, {2, 1}}, {interval, point, interval}), Geometry({0, 2, 5}, 1)};
  const DGConnectivity connectivity(mesh);

  EXPECT_TRUE(connectivity.faces(1).empty());
  for (const std::size_t cell : {0, 2}) {
    const auto faces = connectivity.faces(cell);
    ASSERT_EQ(faces.size(), 2);
    EXPECT_TRUE(faces[0].is_boundary());
    EXPECT_EQ(faces[1].neighbor_cell, 2 - cell);
    EXPECT_EQ(faces[1].neighbor_face, 1);
    EXPECT_EQ(faces[1].vertex_permutation, (std::vector<std::size_t>{0}));
  }
}

TEST(test_dg_connectivity, rejects_a_zero_dimensional_mesh) {
  Mesh mesh{Topology({{0}}, {get_cell_type(CellKind::Point)}), Geometry({0}, 1)};

  EXPECT_THROW(DGConnectivity{mesh}, std::invalid_argument);
}

TEST(test_dg_connectivity, rejects_an_edge_shared_by_three_cells) {
  const auto triangle = get_cell_type(CellKind::Triangle);
  // Three triangles all use the edge between vertices 0 and 1.
  Mesh mesh{Topology({{0, 1, 2}, {1, 0, 3}, {0, 1, 4}}, {triangle, triangle, triangle}),
            Geometry({0, 0, 2, 0, 0, 1, 0, -1, 2, 1}, 2)};

  EXPECT_THROW(DGConnectivity{mesh}, std::invalid_argument);
}

TEST(test_mesh_topology, rejects_two_faces_of_the_same_cell_with_the_same_vertices) {
  const auto interval = get_cell_type(CellKind::Interval);
  Topology topology({{0, 0}}, {interval});

  EXPECT_THROW(topology.calculate_connectivity(), std::invalid_argument);
}

TEST(test_mesh_topology, invalid_face_vertex_preserves_adjacency_from_previous_build) {
  const auto interval = get_cell_type(CellKind::Interval);
  Topology topology({{0, 1}, {1, 2}}, {interval, interval});
  topology.calculate_connectivity();

  topology.conn()[1] = {2};
  EXPECT_THROW(topology.calculate_connectivity(), std::invalid_argument);

  ASSERT_EQ(topology.e_to_e().size(), 2);
  ASSERT_EQ(topology.e_to_f().size(), 2);
  EXPECT_EQ(topology.e_to_e()[0], (std::vector<std::size_t>{0, 1}));
  EXPECT_EQ(topology.e_to_e()[1], (std::vector<std::size_t>{0, 1}));
  EXPECT_EQ(topology.e_to_f()[0], (std::vector<std::size_t>{0, 0}));
  EXPECT_EQ(topology.e_to_f()[1], (std::vector<std::size_t>{1, 1}));
}

TEST(test_dg_connectivity, tetrahedra_store_reciprocal_face_vertex_permutations) {
  const auto tetrahedron = get_cell_type(CellKind::Tetrahedron);
  const auto triangle = get_cell_type(CellKind::Triangle);
  // The tetrahedra share vertices 0, 1, and 2 on z = 0. Cell 1 is a boundary tag.
  Mesh mesh{Topology({{0, 1, 2, 3}, {0, 1, 3}, {1, 2, 0, 4}}, {tetrahedron, triangle, tetrahedron}),
            Geometry({0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1}, 3)};
  const DGConnectivity connectivity(mesh);

  const auto first = connectivity.faces(0);
  const auto second = connectivity.faces(2);
  ASSERT_EQ(first.size(), 4);
  ASSERT_EQ(second.size(), 4);
  EXPECT_TRUE(connectivity.faces(1).empty());
  EXPECT_TRUE(first[0].is_boundary());
  EXPECT_EQ(first[3].neighbor_cell, 2);
  EXPECT_EQ(first[3].neighbor_face, 3);
  EXPECT_EQ(first[3].vertex_permutation, (std::vector<std::size_t>{2, 0, 1}));
  EXPECT_EQ(second[3].neighbor_cell, 0);
  EXPECT_EQ(second[3].neighbor_face, 3);
  EXPECT_EQ(second[3].vertex_permutation, (std::vector<std::size_t>{1, 2, 0}));
}

TEST(test_dg_connectivity, hexahedra_store_reciprocal_face_vertex_permutations) {
  const auto hexahedron = get_cell_type(CellKind::Hexahedron);
  // The first cube's top and the second cube's bottom share vertices 4 through 7.
  Mesh mesh{
      Topology({{0, 1, 2, 3, 4, 5, 6, 7}, {5, 6, 7, 4, 9, 10, 11, 8}}, {hexahedron, hexahedron}),
      Geometry({0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1,
                1, 1, 1, 0, 1, 1, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1, 2},
               3)};
  const DGConnectivity connectivity(mesh);

  const auto first = connectivity.faces(0);
  const auto second = connectivity.faces(1);
  ASSERT_EQ(first.size(), 6);
  ASSERT_EQ(second.size(), 6);
  EXPECT_TRUE(first[0].is_boundary());
  EXPECT_EQ(first[1].neighbor_cell, 1);
  EXPECT_EQ(first[1].neighbor_face, 0);
  EXPECT_EQ(first[1].vertex_permutation, (std::vector<std::size_t>{3, 0, 1, 2}));
  EXPECT_EQ(second[0].neighbor_cell, 0);
  EXPECT_EQ(second[0].neighbor_face, 1);
  EXPECT_EQ(second[0].vertex_permutation, (std::vector<std::size_t>{1, 2, 3, 0}));
}
