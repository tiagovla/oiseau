// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include "oiseau/dg/field.hpp"

using oiseau::dg::Field;

TEST(test_dg_field, retains_its_name) {
  const Field field{"Ez"};

  EXPECT_EQ(field.name(), "Ez");
}

TEST(test_dg_field, copies_share_identity) {
  const Field field{"Ez"};
  const Field copy = field;

  EXPECT_TRUE(field.same_as(copy));
  EXPECT_TRUE(copy.same_as(field));
  EXPECT_EQ(copy.name(), field.name());
}

TEST(test_dg_field, same_named_fields_have_distinct_identity) {
  const Field field{"Ez"};
  const Field another_field{"Ez"};

  EXPECT_EQ(field.name(), another_field.name());
  EXPECT_FALSE(field.same_as(another_field));
}

TEST(test_dg_field, rejects_an_empty_name) {
  EXPECT_THROW((Field{std::string{}}), std::invalid_argument);
}
