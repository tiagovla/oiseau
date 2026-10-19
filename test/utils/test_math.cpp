// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "oiseau/utils/math.hpp"

#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, thresh)                       \
  EXPECT_EQ(expected.size(), actual.size()) << "Array sizes differ.";           \
  for (size_t idx = 0; idx < std::min(expected.size(), actual.size()); ++idx) { \
    EXPECT_NEAR(expected[idx], actual[idx], thresh) << "at index: " << idx;     \
  }

TEST(test_utils, test_jacobi_p) {
  unsigned n = 2;
  double alpha = 1.0, beta = 2.0;
  std::vector<double> input = {0.1, 0.3, 0.7, 0.9};
  std::vector<double> expected = {-0.77365811, -0.6641136, 0.70519279, 1.96495468};
  auto output = oiseau::utils::jacobi_p(n, alpha, beta, input);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}

TEST(test_utils, test_grad_jacobi_p) {
  unsigned n = 2;
  double alpha = 1.0, beta = 2.0;
  std::vector<double> input = {0.1, 0.3, 0.7, 0.9};
  std::vector<double> expected = {-0.41079192, 1.50623703, 5.34029494, 7.25732389};
  auto output = oiseau::utils::grad_jacobi_p(n, alpha, beta, input);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}
