// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "oiseau/io/gmsh_file.hpp"

TEST(test_io, gmsh_parser_mesh_format_handler) {
  std::string str =
      R"(4.1 0 8
$EndMeshFormat)";
  std::stringstream test_stream(str);
  auto s = oiseau::io::detail::mesh_format_handler(test_stream);
  EXPECT_EQ(s.version, 4.1);
  EXPECT_EQ(s.data_size, 8);
  EXPECT_EQ(s.is_binary, 0);
}

TEST(test_io, gmsh_parser_physical_names_handler) {
  std::string str = R"(2
1 11 "region_1"
2 21 "region_2"
$EndPhysicalNames
)";
  std::stringstream test_stream(str);
  auto s = oiseau::io::detail::physical_names_handler(test_stream);
  EXPECT_EQ(s.num_physical_names, 2);
  EXPECT_EQ(s.dimensions[0], 1);
  EXPECT_EQ(s.physical_tags[0], 11);
  EXPECT_EQ(s.names[0], "\"region_1\"");
  EXPECT_EQ(s.num_physical_names, 2);
  EXPECT_EQ(s.dimensions[1], 2);
  EXPECT_EQ(s.physical_tags[1], 21);
  EXPECT_EQ(s.names[1], "\"region_2\"");
}

TEST(test_io, gmsh_parser_entities_handler) {
  std::string str =
      R"(3 0 0 0
1 0 -1 0 0
2 1 -1 0 0
3 1 1 0 0
)";
  std::stringstream test_stream(str);
  auto s = oiseau::io::detail::entities_handler(test_stream, false);
  EXPECT_EQ(s.blocks[0][0].tag, 1);
}

TEST(test_io, gmsh_parser_nodes_handler) {
  std::string str =
      R"(2 2 1 2
0 1 0 1
1
0 -1 0
0 2 0 1
2
1 -1 0)";
  std::stringstream test_stream(str);
  auto s = oiseau::io::detail::nodes_handler(test_stream, false);
  EXPECT_EQ(s.num_entity_blocks, 2);
  EXPECT_EQ(s.max_node_tag, 2);
  EXPECT_EQ(s.min_node_tag, 1);
  EXPECT_EQ(s.num_nodes, 2);
  EXPECT_EQ(s.blocks[0].entity_dim, 0);
  EXPECT_EQ(s.blocks[0].entity_tag, 1);
  EXPECT_EQ(s.blocks[0].parametric, 0);
  EXPECT_EQ(s.blocks[0].node_tags[0], 1);
  EXPECT_EQ(s.blocks[0].node_coords[0], 0);
  EXPECT_EQ(s.blocks[0].node_coords[1], -1);
  EXPECT_EQ(s.blocks[0].node_coords[2], 0);
}
