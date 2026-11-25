// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <stdexcept>

#include "oiseau/mesh/geometry.hpp"

using oiseau::mesh::Geometry;

TEST(test_mesh_geometry, infers_dimension_and_stores_points_in_order) {
  Geometry plane{{1, 2}, {3, 4}};
  EXPECT_EQ(plane.dim(), 2u);
  EXPECT_TRUE(std::ranges::equal(plane.x(), std::array{1.0, 2.0, 3.0, 4.0}));
  EXPECT_TRUE(std::ranges::equal(plane.x_at(1), std::array{3.0, 4.0}));

  const Geometry line{{1}, {2}};
  EXPECT_EQ(line.dim(), 1u);
  EXPECT_TRUE(std::ranges::equal(line.x(), std::array{1.0, 2.0}));
}

TEST(test_mesh_geometry, rejects_empty_and_ragged_point_lists) {
  EXPECT_THROW(Geometry({}), std::invalid_argument);
  EXPECT_THROW(Geometry({{}}), std::invalid_argument);
  EXPECT_THROW(Geometry({{1, 2}, {3}}), std::invalid_argument);
}
