// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#define EXPECT_FLOATS_NEARLY_EQ(expected, actual, threshold)          \
  ASSERT_EQ(expected.size(), actual.size()) << "Array sizes differ."; \
  for (size_t idx = 0; idx < expected.size(); ++idx) {                \
    auto diff = std::abs(expected[idx] - actual[idx]);                \
    if (diff > threshold) {                                           \
      FAIL() << "Stopped at first mismatch:\n"                        \
             << "  index     : " << idx << "\n"                       \
             << "  expected  : " << expected[idx] << "\n"             \
             << "  actual    : " << actual[idx] << "\n"               \
             << "  difference: " << diff << "\n"                      \
             << "  threshold : " << threshold;                        \
    }                                                                 \
  }
