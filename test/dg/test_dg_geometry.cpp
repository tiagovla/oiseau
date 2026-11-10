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
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <xtensor/generators/xbuilder.hpp>

#include "oiseau/dg/detail/geometry.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"

namespace {

oiseau::dg::detail::Point make_point(std::initializer_list<double> values) {
  oiseau::dg::detail::Point point = xt::zeros<double>({values.size()});
  std::ranges::copy(values, point.begin());
  return point;
}

}  // namespace

TEST(test_dg_geometry, vector_norm) {
  EXPECT_DOUBLE_EQ(oiseau::dg::detail::vector_norm(make_point({3.0, 4.0})), 5.0);
}

TEST(test_dg_geometry, padded_matrix) {
  const oiseau::dg::detail::Matrix matrix = {{2.0, 1.0}, {1.0, 3.0}};
  const std::array<double, 9> expected = {2.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0};
  EXPECT_EQ(oiseau::dg::detail::padded_matrix(matrix, 2), expected);
}

TEST(test_dg_geometry, point_array) {
  EXPECT_EQ(oiseau::dg::detail::point_array(make_point({3.0, 4.0})),
            (std::array<double, 3>{3.0, 4.0, 0.0}));
}

TEST(test_dg_geometry, ref_type_for) {
  using oiseau::dg::detail::ref_type_for;
  using oiseau::dg::nodal::RefElementType;
  using oiseau::mesh::CellKind;
  EXPECT_EQ(ref_type_for(CellKind::Interval), RefElementType::Line);
  EXPECT_EQ(ref_type_for(CellKind::Triangle), RefElementType::Triangle);
  EXPECT_EQ(ref_type_for(CellKind::Quadrilateral), RefElementType::Quadrilateral);
  EXPECT_EQ(ref_type_for(CellKind::Tetrahedron), RefElementType::Tetrahedron);
  EXPECT_EQ(ref_type_for(CellKind::Hexahedron), RefElementType::Hexahedron);
  EXPECT_THROW(ref_type_for(CellKind::Undefined), std::invalid_argument);
}

TEST(test_dg_geometry, matrix_scale) {
  const oiseau::dg::detail::Matrix matrix = {{2.0, 1.0}, {1.0, 3.0}};
  EXPECT_DOUBLE_EQ(oiseau::dg::detail::matrix_scale(matrix, 2), 3.0);
}

TEST(test_dg_geometry, is_degenerate) {
  using oiseau::dg::detail::is_degenerate;
  const oiseau::dg::detail::Matrix matrix = {{2.0, 1.0}, {1.0, 3.0}};
  const oiseau::dg::detail::Matrix singular = {{0.0, 0.0}, {0.0, 1.0}};
  EXPECT_FALSE(is_degenerate(matrix, 5.0, 2));
  EXPECT_TRUE(is_degenerate(singular, 0.0, 2));
  EXPECT_TRUE(is_degenerate(matrix, std::numeric_limits<double>::infinity(), 2));
}

TEST(test_dg_geometry, affine_tolerance) {
  const std::vector<oiseau::dg::detail::Point> points = {
      make_point({0.0, 0.0}), make_point({1.0, 0.0}), make_point({0.0, 1.0})};
  const auto tolerance = oiseau::dg::detail::affine_tolerance(points, 2);
  EXPECT_GT(tolerance, 0.0);
  EXPECT_TRUE(std::isfinite(tolerance));
}

TEST(test_dg_geometry, physical_vertex) {
  using oiseau::dg::detail::physical_vertex;
  const oiseau::mesh::Geometry geometry(std::vector<double>{1.0, 2.0, 0.0, 3.0, 4.0, 0.0}, 3);
  const auto vertex = physical_vertex(geometry, 1, 2);
  EXPECT_DOUBLE_EQ(vertex(0), 3.0);
  EXPECT_DOUBLE_EQ(vertex(1), 4.0);
  EXPECT_THROW(physical_vertex(geometry, 2, 2), std::out_of_range);
  const oiseau::mesh::Geometry off_plane(std::vector<double>{0.0, 0.0, 1.0}, 3);
  EXPECT_THROW(physical_vertex(off_plane, 0, 2), std::invalid_argument);
  const oiseau::mesh::Geometry nonfinite(
      std::vector<double>{std::numeric_limits<double>::infinity(), 0.0}, 2);
  EXPECT_THROW(physical_vertex(nonfinite, 0, 2), std::invalid_argument);
}

TEST(test_dg_geometry, independent_vertices) {
  using namespace oiseau::dg::detail;
  const std::vector<Point> triangle = {make_point({0.0, 0.0}), make_point({1.0, 0.0}),
                                       make_point({0.0, 1.0})};
  const auto basis = independent_vertices(triangle, 2);
  EXPECT_EQ(basis[0], 0);
  const auto first = triangle[basis[1]] - triangle[basis[0]];
  const auto second = triangle[basis[2]] - triangle[basis[0]];
  EXPECT_GT(std::abs(first(0) * second(1) - first(1) * second(0)),
            std::numeric_limits<double>::epsilon());
  EXPECT_THROW(independent_vertices(
                   {make_point({0.0, 0.0}), make_point({1.0, 0.0}), make_point({2.0, 0.0})}, 2),
               std::logic_error);
}

TEST(test_dg_geometry, face_measure) {
  using oiseau::dg::detail::face_measure;
  EXPECT_DOUBLE_EQ(face_measure({make_point({2.0})}, 1), 1.0);
  EXPECT_DOUBLE_EQ(face_measure({make_point({0.0, 0.0}), make_point({3.0, 4.0})}, 2), 5.0);
  EXPECT_DOUBLE_EQ(face_measure({make_point({0.0, 0.0, 0.0}), make_point({2.0, 0.0, 0.0}),
                                 make_point({0.0, 3.0, 0.0})},
                                3),
                   3.0);
  EXPECT_DOUBLE_EQ(face_measure({make_point({0.0, 0.0, 0.0}), make_point({1.0, 0.0, 0.0}),
                                 make_point({1.0, 1.0, 0.0}), make_point({0.0, 1.0, 0.0})},
                                3),
                   1.0);
}

TEST(test_dg_geometry, centroid) {
  const std::vector<oiseau::dg::detail::Point> triangle = {
      make_point({0.0, 0.0, 0.0}), make_point({2.0, 0.0, 0.0}), make_point({0.0, 3.0, 0.0})};
  const auto center = oiseau::dg::detail::centroid(triangle);
  EXPECT_DOUBLE_EQ(center(0), 2.0 / 3.0);
  EXPECT_DOUBLE_EQ(center(1), 1.0);
  EXPECT_DOUBLE_EQ(center(2), 0.0);
}

TEST(test_dg_geometry, face_normal) {
  using oiseau::dg::detail::face_normal;
  const auto normal_1d = face_normal({make_point({1.0})}, make_point({0.0}), 1);
  EXPECT_DOUBLE_EQ(normal_1d(0), 1.0);
  const auto normal_2d =
      face_normal({make_point({0.0, 0.0}), make_point({2.0, 0.0})}, make_point({1.0, -1.0}), 2);
  EXPECT_DOUBLE_EQ(normal_2d(0), 0.0);
  EXPECT_DOUBLE_EQ(normal_2d(1), 1.0);
  const auto normal_3d = face_normal(
      {make_point({0.0, 0.0, 0.0}), make_point({2.0, 0.0, 0.0}), make_point({0.0, 3.0, 0.0})},
      make_point({0.0, 0.0, -1.0}), 3);
  EXPECT_DOUBLE_EQ(normal_3d(0), 0.0);
  EXPECT_DOUBLE_EQ(normal_3d(1), 0.0);
  EXPECT_DOUBLE_EQ(normal_3d(2), 1.0);
  EXPECT_THROW(
      face_normal({make_point({0.0, 0.0}), make_point({0.0, 0.0})}, make_point({0.0, -1.0}), 2),
      std::invalid_argument);
}

TEST(test_dg_geometry, intervals_have_signed_jacobians_and_outward_endpoint_normals) {
  const std::string mesh_text = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
1 0 0 3
1
2
3
2 0 0
6 0 0
10 0 0
$EndNodes
$Elements
1 2 1 2
1 0 1 2
1 1 2
2 3 2
$EndElements)";
  const auto mesh = oiseau::io::gmsh_read_from_string(mesh_text);
  const oiseau::dg::DGGeometry geometry(mesh);

  for (std::size_t cell_index = 0; cell_index < 2; ++cell_index) {
    SCOPED_TRACE(cell_index);
    ASSERT_TRUE(geometry.has_cell(cell_index));
    const auto& cell = geometry.cell(cell_index);
    const auto orientation = cell_index == 0 ? 1.0 : -1.0;
    EXPECT_EQ(cell.dimension, 1);
    EXPECT_NEAR(cell.jacobian[0], 2.0 * orientation, 1e-12);
    EXPECT_NEAR(cell.inverse_jacobian[0], 0.5 * orientation, 1e-12);
    for (std::size_t i = 1; i < cell.jacobian.size(); ++i) {
      EXPECT_DOUBLE_EQ(cell.jacobian[i], 0.0);
      EXPECT_DOUBLE_EQ(cell.inverse_jacobian[i], 0.0);
    }
    EXPECT_NEAR(cell.determinant, 2.0 * orientation, 1e-12);
    ASSERT_EQ(cell.faces.size(), 2);
    for (std::size_t face_index = 0; face_index < 2; ++face_index) {
      const auto& face = cell.faces[face_index];
      EXPECT_NEAR(face.outward_normal[0], (face_index == 0 ? -1.0 : 1.0) * orientation, 1e-12);
      EXPECT_DOUBLE_EQ(face.outward_normal[1], 0.0);
      EXPECT_DOUBLE_EQ(face.outward_normal[2], 0.0);
      EXPECT_DOUBLE_EQ(face.measure, 1.0);
      EXPECT_DOUBLE_EQ(face.surface_jacobian, 1.0);
    }
  }
}

TEST(test_dg_geometry, skew_triangle_has_padded_jacobians_and_outward_normals) {
  const std::string mesh_text = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
2 0 0 3
1
2
3
0 0 0
2 1 0
1 3 0
$EndNodes
$Elements
1 1 1 1
2 0 2 1
1 1 2 3
$EndElements)";
  const auto mesh = oiseau::io::gmsh_read_from_string(mesh_text);
  const oiseau::dg::DGGeometry geometry(mesh);

  ASSERT_TRUE(geometry.has_cell(0));
  const auto& cell = geometry.cell(0);
  EXPECT_EQ(cell.dimension, 2);
  const std::array<double, 9> expected_jacobian = {1.0, 0.5, 0.0, 0.5, 1.5, 0.0, 0.0, 0.0, 0.0};
  const std::array<double, 9> expected_inverse = {1.2, -0.4, 0.0, -0.4, 0.8, 0.0, 0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < 9; ++i) {
    EXPECT_NEAR(cell.jacobian[i], expected_jacobian[i], 1e-12);
    EXPECT_NEAR(cell.inverse_jacobian[i], expected_inverse[i], 1e-12);
  }
  EXPECT_NEAR(cell.determinant, 1.25, 1e-12);

  const double inverse_sqrt_five = 1.0 / std::sqrt(5.0);
  const double inverse_sqrt_ten = 1.0 / std::sqrt(10.0);
  const std::array<std::array<double, 3>, 3> expected_normals = {
      std::array<double, 3>{2.0 * inverse_sqrt_five, inverse_sqrt_five, 0.0},
      std::array<double, 3>{-3.0 * inverse_sqrt_ten, inverse_sqrt_ten, 0.0},
      std::array<double, 3>{inverse_sqrt_five, -2.0 * inverse_sqrt_five, 0.0},
  };
  ASSERT_EQ(cell.faces.size(), expected_normals.size());
  for (std::size_t face = 0; face < cell.faces.size(); ++face) {
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      EXPECT_NEAR(cell.faces[face].outward_normal[coordinate], expected_normals[face][coordinate],
                  1e-12);
    }
  }
}

TEST(test_dg_geometry, one_tetrahedron_has_outward_face_normals) {
  const std::string mesh_text = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 0 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 0 4 1
1 1 2 3 4
$EndElements)";
  const auto mesh = oiseau::io::gmsh_read_from_string(mesh_text);
  const oiseau::dg::DGGeometry geometry(mesh);

  ASSERT_TRUE(geometry.has_cell(0));
  const auto& cell = geometry.cell(0);
  ASSERT_EQ(cell.faces.size(), 4);

  const double inverse_sqrt_three = 1.0 / std::sqrt(3.0);
  const std::array<std::array<double, 3>, 4> expected_normals = {
      std::array<double, 3>{inverse_sqrt_three, inverse_sqrt_three, inverse_sqrt_three},
      std::array<double, 3>{-1.0, 0.0, 0.0},
      std::array<double, 3>{0.0, -1.0, 0.0},
      std::array<double, 3>{0.0, 0.0, -1.0},
  };
  const std::array<double, 4> expected_measures = {
      std::sqrt(3.0) / 2.0,
      0.5,
      0.5,
      0.5,
  };

  for (std::size_t face_index = 0; face_index < cell.faces.size(); ++face_index) {
    const auto& face = cell.faces[face_index];
    EXPECT_NEAR(face.measure, expected_measures[face_index], 1e-12);
    EXPECT_GT(face.surface_jacobian, 0.0);
    for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
      EXPECT_NEAR(face.outward_normal[coordinate], expected_normals[face_index][coordinate], 1e-12);
    }
  }
}
