// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>

#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/io/gmsh.hpp"

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
