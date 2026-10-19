// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

TEST(test_io, gmsh_read_from_string_3d_tetra_block) {
  std::string str =
      R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 8 1 8
3 0 0 8
1
2
3
4
5
6
7
8
0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00
1.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00
1.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00
0.0000000000000000e+00 0.0000000000000000e+00 1.0000000000000000e+00
1.0000000000000000e+00 0.0000000000000000e+00 1.0000000000000000e+00
1.0000000000000000e+00 1.0000000000000000e+00 1.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+00 1.0000000000000000e+00
$EndNodes
$Elements
1 5 1 5
3 0 4 5
1 1 2 4 5
2 2 3 4 7
3 2 4 5 7
4 2 5 6 7
5 4 5 7 8
$EndElements)";
  oiseau::mesh::Mesh mesh = oiseau::io::gmsh_read_from_string(str);
  auto conn = mesh.topology().conn();
  std::vector<std::vector<size_t>> actual(conn.begin(), conn.end());
  std::vector<std::vector<size_t>> expected = {
      {0, 1, 3, 4}, {1, 2, 3, 6}, {1, 3, 4, 6}, {1, 4, 5, 6}, {3, 4, 6, 7}};
  EXPECT_EQ(actual, expected);
}

TEST(test_io, gmsh_read_from_string_3d_tetra_hexa_block) {
  std::string str =
      R"($MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 2
1 0 0 0 0 0 0 1 1 0
2 0 0 0 0 0 0 1 2 0
$EndEntities
$Nodes
2 12 1 12
3 1 0 8
1
2
3
4
5
6
7
8
0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00
5.0000000000000000e-01 0.0000000000000000e+00 0.0000000000000000e+00
5.0000000000000000e-01 1.0000000000000000e+00 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00
0.0000000000000000e+00 0.0000000000000000e+00 1.0000000000000000e+00
5.0000000000000000e-01 0.0000000000000000e+00 1.0000000000000000e+00
5.0000000000000000e-01 1.0000000000000000e+00 1.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+00 1.0000000000000000e+00
3 2 0 4
9
10
11
12
1.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00
1.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00
1.0000000000000000e+00 0.0000000000000000e+00 1.0000000000000000e+00
1.0000000000000000e+00 1.0000000000000000e+00 1.0000000000000000e+00
$EndNodes
$Elements
2 3 1 3
3 1 5 1
1 1 2 3 4 5 6 7 8
3 2 4 2
2 2 9 3 6
3 6 3 9 12
$EndElements)";
  oiseau::mesh::Mesh mesh = oiseau::io::gmsh_read_from_string(str);
  auto conn = mesh.topology().conn();
  std::vector<std::vector<size_t>> actual(conn.begin(), conn.end());
  std::vector<std::vector<size_t>> expected = {
      {0, 1, 2, 3, 4, 5, 6, 7}, {1, 8, 2, 5}, {5, 2, 8, 11}};
  EXPECT_EQ(actual, expected);
}

TEST(test_io, gmsh_celltype_to_oiseau_celltype) {
  EXPECT_EQ(oiseau::io::detail::gmsh_celltype_to_oiseau_celltype(15),
            oiseau::mesh::get_cell_type(oiseau::mesh::CellKind::Point));
  EXPECT_EQ(oiseau::io::detail::gmsh_celltype_to_oiseau_celltype(2),
            oiseau::mesh::get_cell_type(oiseau::mesh::CellKind::Triangle));
  EXPECT_THROW(oiseau::io::detail::gmsh_celltype_to_oiseau_celltype(420), std::runtime_error);
}
