// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <utility>
#include <vector>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/core/xshape.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/utils/mdarray.hpp"

namespace oiseau {
namespace {
TEST(test_mdspan, xtensor_move_to_mdspan) {
  using std_xtensor = xt::xtensor_container<std::vector<double>, 2>;
  using exts_type = std::extents<std::size_t, 5, 2>;
  std_xtensor tensor = std_xtensor::from_shape({5, 2});
  double *ptr_tensor = tensor.data();
  exts_type exts;
  stdex::mdarray<double, exts_type> mdarray(exts, std::move(tensor.storage()));
  double *ptr_mdarray = mdarray.data();
  EXPECT_EQ(ptr_tensor, ptr_mdarray);
}
}  // namespace
}  // namespace oiseau
