// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/region.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau::mesh {
namespace {

TEST(MeshRegions, MergesRepeatedNamesAndAllowsOverlap) {
  Mesh mesh{Topology({{0, 1}, {1, 2}, {2, 3}},
                     {CellKind::Interval, CellKind::Interval, CellKind::Interval}),
            Geometry({{0}, {1}, {2}, {3}})};

  mesh.add_region("dielectric", {2, 1});
  const auto dielectric_tag = mesh.region("dielectric").tag();
  mesh.add_region("dielectric", {1, 0});
  mesh.add_region("pml", {2});

  EXPECT_EQ(mesh.region("dielectric").name(), "dielectric");
  EXPECT_EQ(mesh.region("dielectric").tag(), dielectric_tag);
  EXPECT_NE(mesh.region("pml").tag(), dielectric_tag);
  EXPECT_EQ(mesh.region("dielectric").cells(), (std::vector<std::size_t>{0, 1, 2}));
  EXPECT_EQ(mesh.region("pml").cells(), (std::vector<std::size_t>{2}));
}

TEST(MeshRegions, RejectsInvalidAssignmentsBeforeChangingMembership) {
  Mesh mesh{
      Topology({{0, 1}, {1, 2}, {2}}, {CellKind::Interval, CellKind::Interval, CellKind::Point}),
      Geometry({{0}, {1}, {2}})};
  mesh.add_region("dielectric", {0});

  EXPECT_THROW(mesh.add_region("dielectric", {1, 3}), std::out_of_range);
  EXPECT_THROW(mesh.add_region("dielectric", {1, 2}), std::invalid_argument);
  EXPECT_THROW(mesh.add_region("", {1}), std::invalid_argument);
  EXPECT_THROW(mesh.region("pml"), std::out_of_range);
  EXPECT_EQ(mesh.region("dielectric").cells(), (std::vector<std::size_t>{0}));
}

TEST(Region, UnionRequiresMatchingTagAndName) {
  const Region left{7, "dielectric", {3, 1}};
  const Region right{7, "dielectric", {4, 1}};
  const Region other_tag{8, "dielectric", {4}};
  const Region other_name{7, "pml", {4}};

  const auto merged = left | right;

  EXPECT_EQ(merged.tag(), 7u);
  EXPECT_EQ(merged.name(), "dielectric");
  EXPECT_EQ(merged.cells(), (std::vector<std::size_t>{1, 3, 4}));
  EXPECT_EQ(left.cells(), (std::vector<std::size_t>{1, 3}));
  EXPECT_EQ(right.cells(), (std::vector<std::size_t>{1, 4}));
  EXPECT_THROW(left | other_tag, std::invalid_argument);
  EXPECT_THROW(left | other_name, std::invalid_argument);
}

}  // namespace
}  // namespace oiseau::mesh
