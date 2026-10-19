// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cstddef>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oiseau/utils/jagged_array.hpp"

using oiseau::utils::JaggedArray;

#include "gtest/gtest.h"

class jagged_array_test_fixture
    : public ::testing::Test {  // Fixture name also in snake_case for consistency
 protected:
  JaggedArray<int> arr_int_empty;
  JaggedArray<std::size_t> arr_size_t_empty;

  JaggedArray<int> arr_int_populated = {{1, 2}, {3, 4, 5}, {6}};
  JaggedArray<std::size_t> arr_size_t_populated = {{0, 1}, {10}, {20, 21, 22}};
};

TEST(jagged_array_constructor, default_constructor) {
  JaggedArray<int> ja;
  EXPECT_TRUE(ja.is_empty());
  EXPECT_EQ(ja.num_rows(), 0);
  EXPECT_EQ(ja.total_elements(), 0);
}

TEST(jagged_array_constructor, num_rows_constructor) {
  JaggedArray<int> ja(3);
  EXPECT_FALSE(ja.is_empty());
  EXPECT_EQ(ja.num_rows(), 3);
  EXPECT_EQ(ja.total_elements(), 0);
  for (std::size_t i = 0; i < 3; ++i) {
    EXPECT_EQ(ja.num_cols(i), 0);
  }
}

TEST(jagged_array_constructor, initializer_list_constructor_int) {
  JaggedArray<int> ja = {{1, 2, 3}, {4, 5}, {}, {6, 7, 8, 9}};
  EXPECT_FALSE(ja.is_empty());
  EXPECT_EQ(ja.num_rows(), 4);
  EXPECT_EQ(ja.total_elements(), 9);

  EXPECT_EQ(ja.num_cols(0), 3);
  EXPECT_EQ(ja.at(0, 0), 1);
  EXPECT_EQ(ja.at(0, 1), 2);
  EXPECT_EQ(ja.at(0, 2), 3);

  EXPECT_EQ(ja.num_cols(1), 2);
  EXPECT_EQ(ja.at(1, 0), 4);
  EXPECT_EQ(ja.at(1, 1), 5);

  EXPECT_EQ(ja.num_cols(2), 0);

  EXPECT_EQ(ja.num_cols(3), 4);
  EXPECT_EQ(ja.at(3, 0), 6);
  EXPECT_EQ(ja.at(3, 3), 9);
}

TEST(jagged_array_constructor, initializer_list_constructor_size_t) {
  JaggedArray<std::size_t> ja = {{100, 200}, {300}};
  EXPECT_EQ(ja.num_rows(), 2);
  EXPECT_EQ(ja.num_cols(0), 2);
  EXPECT_EQ(ja.num_cols(1), 1);
  EXPECT_EQ(ja.at(0, 0), 100u);
  EXPECT_EQ(ja.at(1, 0), 300u);
}

TEST(jagged_array_constructor, initializer_list_constructor_vector_t) {
  std::vector<int> r1 = {10, 20};
  std::vector<int> r2 = {};
  std::vector<int> r3 = {30, 40, 50};
  JaggedArray<int> ja = {r1, r2, r3};

  EXPECT_EQ(ja.num_rows(), 3);
  EXPECT_EQ(ja.num_cols(0), 2);
  EXPECT_EQ(ja.at(0, 1), 20);
  EXPECT_EQ(ja.num_cols(1), 0);
  EXPECT_EQ(ja.num_cols(2), 3);
  EXPECT_EQ(ja.at(2, 2), 50);
}

TEST(jagged_array_constructor, copy_constructor) {
  JaggedArray<int> original = {{1, 2}, {3}};
  JaggedArray<int> copy = original;

  EXPECT_EQ(copy.num_rows(), 2);
  EXPECT_EQ(copy.num_cols(0), 2);
  EXPECT_EQ(copy.at(0, 0), 1);
  EXPECT_EQ(copy.at(1, 0), 3);

  copy.add_element(0, 99);
  EXPECT_EQ(original.at(0, 1), 2);
  EXPECT_EQ(copy.at(0, 2), 99);
}

TEST(jagged_array_constructor, move_constructor) {
  JaggedArray<int> original = {{1, 2}, {3}};
  JaggedArray<int> moved_to = std::move(original);

  EXPECT_EQ(moved_to.num_rows(), 2);
  EXPECT_EQ(moved_to.at(0, 0), 1);

  EXPECT_TRUE(original.is_empty() || original.num_rows() == 0);
}

TEST(jagged_array_assignment, copy_assignment) {
  JaggedArray<int> ja1 = {{1}, {2, 3}};
  JaggedArray<int> ja2 = {{4, 5, 6}};
  ja2 = ja1;

  EXPECT_EQ(ja2.num_rows(), 2);
  EXPECT_EQ(ja2.at(1, 1), 3);

  ja1.add_element(0, 9);
  EXPECT_EQ(ja2.at(0, 0), 1);

  ja1 = ja1;
  EXPECT_EQ(ja1.num_rows(), 2);
  EXPECT_EQ(ja1.at(0, 1), 9);
}

TEST(jagged_array_assignment, move_assignment) {
  JaggedArray<int> ja1 = {{1}, {2, 3}};
  JaggedArray<int> ja2 = {{4, 5, 6}};
  ja2 = std::move(ja1);

  EXPECT_EQ(ja2.num_rows(), 2);
  EXPECT_EQ(ja2.at(1, 1), 3);
  EXPECT_TRUE(ja1.is_empty() || ja1.num_rows() == 0);
}

TEST_F(jagged_array_test_fixture, accessors_num_rows_cols_total) {
  EXPECT_EQ(arr_int_empty.num_rows(), 0);
  EXPECT_EQ(arr_int_empty.total_elements(), 0);

  EXPECT_EQ(arr_int_populated.num_rows(), 3);
  EXPECT_EQ(arr_int_populated.total_elements(), 6);
  EXPECT_EQ(arr_int_populated.num_cols(0), 2);
  EXPECT_EQ(arr_int_populated.num_cols(1), 3);
  EXPECT_EQ(arr_int_populated.num_cols(2), 1);
  EXPECT_THROW(arr_int_populated.num_cols(3), std::out_of_range);
}

TEST_F(jagged_array_test_fixture, accessor_operator_square_brackets) {
  std::span<int> row0 = arr_int_populated[0];
  ASSERT_EQ(row0.size(), 2u);
  EXPECT_EQ(row0[0], 1);
  EXPECT_EQ(row0[1], 2);

  arr_int_populated[0][0] = 100;
  EXPECT_EQ(arr_int_populated.at(0, 0), 100);

  const JaggedArray<int>& const_arr = arr_int_populated;
  std::span<const int> const_row1 = const_arr[1];
  ASSERT_EQ(const_row1.size(), 3u);
  EXPECT_EQ(const_row1[0], 3);

  EXPECT_THROW(arr_int_populated[3], std::out_of_range);
  EXPECT_THROW(const_arr[3], std::out_of_range);
  EXPECT_THROW(arr_int_empty[0], std::out_of_range);
}

TEST_F(jagged_array_test_fixture, accessor_at) {
  EXPECT_EQ(arr_int_populated.at(0, 0), 1);
  EXPECT_EQ(arr_int_populated.at(1, 2), 5);
  arr_int_populated.at(1, 2) = 50;
  EXPECT_EQ(arr_int_populated.at(1, 2), 50);

  const JaggedArray<int>& const_arr = arr_int_populated;
  EXPECT_EQ(const_arr.at(0, 1), 2);

  EXPECT_THROW(arr_int_populated.at(3, 0), std::out_of_range);
  EXPECT_THROW(arr_int_populated.at(0, 2), std::out_of_range);
  EXPECT_THROW(const_arr.at(1, 3), std::out_of_range);
}

TEST(jagged_array_modifier, add_row_initializer_list) {
  JaggedArray<int> ja;
  ja.add_row({10, 20});
  EXPECT_EQ(ja.num_rows(), 1u);
  EXPECT_EQ(ja.num_cols(0), 2u);
  EXPECT_EQ(ja.at(0, 1), 20);

  ja.add_row({30});
  EXPECT_EQ(ja.num_rows(), 2u);
  EXPECT_EQ(ja.num_cols(1), 1u);
  EXPECT_EQ(ja.at(1, 0), 30);
  EXPECT_EQ(ja.total_elements(), 3u);
}

TEST(jagged_array_modifier, add_row_const_vector) {
  JaggedArray<int> ja;
  std::vector<int> row_data = {11, 22};
  ja.add_row(row_data);
  EXPECT_EQ(ja.num_rows(), 1u);
  EXPECT_EQ(ja.at(0, 0), 11);
}

TEST(jagged_array_modifier, add_row_move_vector) {
  JaggedArray<int> ja;
  std::vector<int> row_data = {77, 88};
  ja.add_row(std::move(row_data));
  EXPECT_EQ(ja.num_rows(), 1u);
  EXPECT_EQ(ja.at(0, 0), 77);
  EXPECT_TRUE(row_data.empty());
}

TEST(jagged_array_modifier, insert_row) {
  JaggedArray<int> ja = {{10}, {40, 50}};
  ja.insert_row(1, {20, 30});
  ASSERT_EQ(ja.num_rows(), 3u);
  EXPECT_EQ(ja.at(0, 0), 10);
  EXPECT_EQ(ja.at(1, 0), 20);
  EXPECT_EQ(ja.at(1, 1), 30);
  EXPECT_EQ(ja.at(2, 0), 40);
  EXPECT_EQ(ja.total_elements(), 5u);

  ja.insert_row(0, {5});
  ASSERT_EQ(ja.num_rows(), 4u);
  EXPECT_EQ(ja.at(0, 0), 5);
  EXPECT_EQ(ja.at(1, 0), 10);

  ja.insert_row(ja.num_rows(), {60});
  ASSERT_EQ(ja.num_rows(), 5u);
  EXPECT_EQ(ja.at(4, 0), 60);

  EXPECT_THROW(ja.insert_row(10, {0}), std::out_of_range);
}

TEST(jagged_array_modifier, remove_row) {
  JaggedArray<int> ja = {{1}, {2, 3}, {4, 5, 6}, {7}};

  ja.remove_row(1);
  ASSERT_EQ(ja.num_rows(), 3u);
  EXPECT_EQ(ja.at(0, 0), 1);
  EXPECT_EQ(ja.at(1, 0), 4);
  EXPECT_EQ(ja.at(1, 2), 6);
  EXPECT_EQ(ja.at(2, 0), 7);
  EXPECT_EQ(ja.total_elements(), 5u);

  ja.remove_row(ja.num_rows() - 1);
  ASSERT_EQ(ja.num_rows(), 2u);
  EXPECT_EQ(ja.at(1, 0), 4);
  EXPECT_EQ(ja.total_elements(), 4u);

  ja.remove_row(0);
  ASSERT_EQ(ja.num_rows(), 1u);
  EXPECT_EQ(ja.at(0, 0), 4);
  EXPECT_EQ(ja.total_elements(), 3u);

  ja.remove_row(0);
  EXPECT_TRUE(ja.is_empty());
  EXPECT_EQ(ja.total_elements(), 0u);

  EXPECT_THROW(ja.remove_row(0), std::out_of_range);
}

TEST(jagged_array_modifier, add_element) {
  JaggedArray<int> ja = {{1}, {}};
  ja.add_element(0, 2);
  EXPECT_EQ(ja.num_cols(0), 2u);
  EXPECT_EQ(ja.at(0, 1), 2);

  ja.add_element(1, 3);
  EXPECT_EQ(ja.num_cols(1), 1u);
  EXPECT_EQ(ja.at(1, 0), 3);
  EXPECT_EQ(ja.total_elements(), 3u);

  int val = 4;
  ja.add_element(1, std::move(val));
  EXPECT_EQ(ja.num_cols(1), 2u);
  EXPECT_EQ(ja.at(1, 1), 4);

  EXPECT_THROW(ja.add_element(2, 5), std::out_of_range);
}

TEST(jagged_array_modifier, clear) {
  JaggedArray<int> ja = {{1, 2}, {3}};
  ja.clear();
  EXPECT_TRUE(ja.is_empty());
  EXPECT_EQ(ja.num_rows(), 0u);
  EXPECT_EQ(ja.total_elements(), 0u);

  JaggedArray<int> empty_ja;
  empty_ja.clear();
  EXPECT_TRUE(empty_ja.is_empty());
}

TEST_F(jagged_array_test_fixture, is_empty) {
  EXPECT_TRUE(arr_int_empty.is_empty());
  EXPECT_FALSE(arr_int_populated.is_empty());
  arr_int_populated.clear();
  EXPECT_TRUE(arr_int_populated.is_empty());
}

TEST_F(jagged_array_test_fixture, iterator_empty_range) {
  int count = 0;
  for (auto row_span : arr_int_empty) {
    (void)row_span;
    count++;
  }
  EXPECT_EQ(count, 0);

  EXPECT_EQ(arr_int_empty.begin(), arr_int_empty.end());
  EXPECT_EQ(arr_int_empty.cbegin(), arr_int_empty.cend());
}

TEST_F(jagged_array_test_fixture, iterator_read_and_modify) {
  std::vector<std::vector<int>> expected_values = {{1, 2}, {3, 4, 5}, {6}};
  std::size_t row_idx = 0;
  for (auto row_span : arr_int_populated) {
    ASSERT_LT(row_idx, expected_values.size());
    ASSERT_EQ(row_span.size(), expected_values[row_idx].size());
    for (std::size_t col_idx = 0; col_idx < row_span.size(); ++col_idx) {
      EXPECT_EQ(row_span[col_idx], expected_values[row_idx][col_idx]);
      row_span[col_idx] *= 2;
    }
    row_idx++;
  }
  EXPECT_EQ(row_idx, 3u);

  EXPECT_EQ(arr_int_populated.at(0, 0), 2);
  EXPECT_EQ(arr_int_populated.at(1, 1), 8);
  EXPECT_EQ(arr_int_populated.at(2, 0), 12);
}

TEST_F(jagged_array_test_fixture, const_iterator_read) {
  const JaggedArray<std::size_t>& const_arr = arr_size_t_populated;
  std::vector<std::vector<std::size_t>> expected_values = {{0, 1}, {10}, {20, 21, 22}};
  std::size_t row_idx = 0;
  for (std::span<const std::size_t> row_span : const_arr) {
    ASSERT_LT(row_idx, expected_values.size());
    ASSERT_EQ(row_span.size(), expected_values[row_idx].size());
    for (std::size_t col_idx = 0; col_idx < row_span.size(); ++col_idx) {
      EXPECT_EQ(row_span[col_idx], expected_values[row_idx][col_idx]);
    }
    row_idx++;
  }
  EXPECT_EQ(row_idx, 3u);
}

TEST_F(jagged_array_test_fixture, iterator_operations) {
  auto it = arr_int_populated.begin();
  auto it_end = arr_int_populated.end();

  ASSERT_NE(it, it_end);
  EXPECT_EQ((*it).size(), 2u);

  ++it;
  ASSERT_NE(it, it_end);
  EXPECT_EQ((*it).size(), 3u);

  auto it2 = arr_int_populated.begin();
  it2 += 2;
  ASSERT_NE(it2, it_end);
  EXPECT_EQ((*it2).size(), 1u);
  EXPECT_EQ((*it2)[0], 6);

  EXPECT_EQ(it_end - arr_int_populated.begin(), 3);
  EXPECT_EQ((arr_int_populated.begin() + 1), it);

  it = it_end;
  --it;
  EXPECT_EQ((*it)[0], 6);

  auto end_it = arr_int_populated.end();
  std::span<int> end_span = *end_it;
  EXPECT_TRUE(end_span.empty());

  auto empty_arr_end_it = arr_int_empty.end();
  std::span<int> empty_arr_end_span = *empty_arr_end_it;
  EXPECT_TRUE(empty_arr_end_span.empty());
}

TEST(jagged_array_stream_output, output_format) {
  JaggedArray<int> empty_ja;
  std::stringstream ss_empty;
  ss_empty << empty_ja;
  std::string expected_empty = "JaggedArray (0 rows, 0 total elements):\n  (empty)\n";
  EXPECT_EQ(ss_empty.str(), expected_empty);

  JaggedArray<std::size_t> ja = {{1, 2}, {}, {3, 4, 5}};
  std::stringstream ss;
  ss << ja;
  std::string expected_non_empty =
      "JaggedArray (3 rows, 5 total elements):\n"
      "  [0] (2 elements): {1, 2}\n"
      "  [1] (0 elements): {}\n"
      "  [2] (3 elements): {3, 4, 5}\n";
  EXPECT_EQ(ss.str(), expected_non_empty);
}

TEST(jagged_array_integer, complex_operations) {
  JaggedArray<int> ja;
  ja.add_row({10, 20});
  ja.add_row({});
  ja.add_row({30});

  ASSERT_EQ(ja.num_rows(), 3u);
  ja.insert_row(1, {100, 200});
  ASSERT_EQ(ja.num_rows(), 4u);
  EXPECT_EQ(ja.at(0, 1), 20);
  EXPECT_EQ(ja.at(1, 0), 100);
  EXPECT_EQ(ja.num_cols(2), 0u);
  EXPECT_EQ(ja.at(3, 0), 30);

  ja.remove_row(2);
  ASSERT_EQ(ja.num_rows(), 3u);
  EXPECT_EQ(ja.at(2, 0), 30);

  ja.add_element(0, 300);
  EXPECT_EQ(ja[0].size(), 3u);
  EXPECT_EQ(ja.at(0, 2), 300);
  EXPECT_EQ(ja.total_elements(), 6u);
}
