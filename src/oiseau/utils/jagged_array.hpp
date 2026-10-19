// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <algorithm>         // For std::copy, std::move, std::next
#include <initializer_list>  // For std::initializer_list
#include <iostream>          // For std::ostream and operator<<
#include <iterator>          // For std::iterator related tags
#include <span>              // For std::span
#include <stdexcept>         // For std::out_of_range
#include <string>            // For std::to_string in error messages
#include <vector>

namespace oiseau::utils {

template <typename T>
class JaggedArray;

template <typename T>
std::ostream& operator<<(std::ostream& os, const JaggedArray<T>& arr);

template <typename T>
class JaggedArray {
 private:
  std::vector<T> m_data;
  std::vector<std::size_t> m_row_offsets;

  void check_row_bounds(std::size_t r_idx, const char* func_name) const {
    if (r_idx >= num_rows()) {
      throw std::out_of_range(std::string(func_name) + " - Row index (" + std::to_string(r_idx) +
                              ") out of bounds for " + std::to_string(num_rows()) + " rows.");
    }
  }

  void check_row_insert_bounds(std::size_t r_idx, const char* func_name) const {
    if (r_idx > num_rows()) {
      throw std::out_of_range(std::string(func_name) + " - Row index for insertion (" +
                              std::to_string(r_idx) + ") out of bounds for " +
                              std::to_string(num_rows()) + " rows.");
    }
  }

  static std::size_t get_total_elements_hint(
      const std::initializer_list<std::initializer_list<T>>& il) {
    std::size_t total = 0;
    for (const auto& inner : il) total += inner.size();
    return total;
  }
  static std::size_t get_total_elements_hint(const std::initializer_list<std::vector<T>>& il) {
    std::size_t total = 0;
    for (const auto& inner : il) total += inner.size();
    return total;
  }

 public:
  JaggedArray() { m_row_offsets.push_back(0); }

  explicit JaggedArray(std::size_t num_r) : m_row_offsets(num_r + 1, 0) {}

  JaggedArray(std::initializer_list<std::initializer_list<T>> il) {
    m_row_offsets.push_back(0);
    std::size_t current_offset = 0;
    m_data.reserve(get_total_elements_hint(il));
    for (const auto& inner_list : il) {
      m_data.insert(m_data.end(), inner_list.begin(), inner_list.end());
      current_offset += inner_list.size();
      m_row_offsets.push_back(current_offset);
    }
  }

  JaggedArray(std::initializer_list<std::vector<T>> il) {
    m_row_offsets.push_back(0);
    std::size_t current_offset = 0;
    m_data.reserve(get_total_elements_hint(il));
    for (const auto& row_vec : il) {
      m_data.insert(m_data.end(), row_vec.begin(), row_vec.end());
      current_offset += row_vec.size();
      m_row_offsets.push_back(current_offset);
    }
  }

  JaggedArray(const JaggedArray& other) = default;
  JaggedArray(JaggedArray&& other) noexcept = default;
  JaggedArray& operator=(const JaggedArray& other) = default;
  JaggedArray& operator=(JaggedArray&& other) noexcept = default;

  std::size_t num_rows() const noexcept {
    if (m_row_offsets.empty()) return 0;
    return m_row_offsets.size() - 1;
  }

  std::size_t num_cols(std::size_t r_idx) const {
    check_row_bounds(r_idx, "num_cols");
    return m_row_offsets[r_idx + 1] - m_row_offsets[r_idx];
  }

  std::size_t total_elements() const noexcept { return m_data.size(); }

  std::span<T> operator[](std::size_t r_idx) {
    if (r_idx >= num_rows()) {
      throw std::out_of_range("JaggedArray::operator[] - Row index (" + std::to_string(r_idx) +
                              ") out of bounds.");
    }
    std::size_t start_offset = m_row_offsets[r_idx];
    std::size_t length = m_row_offsets[r_idx + 1] - start_offset;
    return std::span<T>(m_data.data() + start_offset, length);
  }

  std::span<const T> operator[](std::size_t r_idx) const {
    if (r_idx >= num_rows()) {
      throw std::out_of_range("JaggedArray::operator[] const - Row index (" +
                              std::to_string(r_idx) + ") out of bounds.");
    }
    std::size_t start_offset = m_row_offsets[r_idx];
    std::size_t length = m_row_offsets[r_idx + 1] - start_offset;
    return std::span<const T>(m_data.data() + start_offset, length);
  }

  T& at(std::size_t r_idx, std::size_t c_idx) {
    check_row_bounds(r_idx, "at (row)");
    std::size_t start_offset = m_row_offsets[r_idx];
    std::size_t row_len = m_row_offsets[r_idx + 1] - start_offset;
    if (c_idx >= row_len) {
      throw std::out_of_range("JaggedArray::at - Column index (" + std::to_string(c_idx) +
                              ") out of bounds for row " + std::to_string(r_idx) + " with length " +
                              std::to_string(row_len) + ".");
    }
    return m_data[start_offset + c_idx];
  }

  const T& at(std::size_t r_idx, std::size_t c_idx) const {
    check_row_bounds(r_idx, "at (row)");
    std::size_t start_offset = m_row_offsets[r_idx];
    std::size_t row_len = m_row_offsets[r_idx + 1] - start_offset;
    if (c_idx >= row_len) {
      throw std::out_of_range("JaggedArray::at - Column index (" + std::to_string(c_idx) +
                              ") out of bounds for row " + std::to_string(r_idx) + " with length " +
                              std::to_string(row_len) + ".");
    }
    return m_data[start_offset + c_idx];
  }

  void add_row(std::initializer_list<T> il) {
    m_data.insert(m_data.end(), il.begin(), il.end());
    m_row_offsets.push_back(m_data.size());
  }

  void add_row(const std::vector<T>& new_row_data) {
    m_data.insert(m_data.end(), new_row_data.begin(), new_row_data.end());
    m_row_offsets.push_back(m_data.size());
  }

  void add_row(std::vector<T>&& new_row_data) {
    m_data.insert(m_data.end(), std::make_move_iterator(new_row_data.begin()),
                  std::make_move_iterator(new_row_data.end()));
    m_row_offsets.push_back(m_data.size());
    new_row_data.clear();
  }

  void insert_row(std::size_t r_idx, std::initializer_list<T> il) {
    check_row_insert_bounds(r_idx, "insert_row");

    std::size_t insert_data_pos = m_row_offsets[r_idx];
    std::size_t new_row_len = il.size();

    m_data.insert(m_data.begin() + insert_data_pos, il.begin(), il.end());
    m_row_offsets.insert(m_row_offsets.begin() + r_idx + 1, m_row_offsets[r_idx] + new_row_len);

    for (std::size_t i = r_idx + 2; i < m_row_offsets.size(); ++i) {
      m_row_offsets[i] += new_row_len;
    }
  }

  void insert_row(std::size_t r_idx, const std::vector<T>& new_row_data) {
    check_row_insert_bounds(r_idx, "insert_row");

    std::size_t insert_data_pos = m_row_offsets[r_idx];
    std::size_t new_row_len = new_row_data.size();

    m_data.insert(m_data.begin() + insert_data_pos, new_row_data.begin(), new_row_data.end());
    m_row_offsets.insert(m_row_offsets.begin() + r_idx + 1, m_row_offsets[r_idx] + new_row_len);
    for (std::size_t i = r_idx + 2; i < m_row_offsets.size(); ++i) {
      m_row_offsets[i] += new_row_len;
    }
  }

  void remove_row(std::size_t r_idx) {
    check_row_bounds(r_idx, "remove_row");

    std::size_t start_data_pos = m_row_offsets[r_idx];
    std::size_t end_data_pos = m_row_offsets[r_idx + 1];
    std::size_t removed_len = end_data_pos - start_data_pos;

    m_data.erase(m_data.begin() + start_data_pos, m_data.begin() + end_data_pos);
    m_row_offsets.erase(m_row_offsets.begin() + r_idx + 1);

    for (std::size_t i = r_idx + 1; i < m_row_offsets.size(); ++i) {
      m_row_offsets[i] -= removed_len;
    }
  }

  void add_element(std::size_t r_idx, const T& value) {
    check_row_bounds(r_idx, "add_element");

    std::size_t insert_data_pos = m_row_offsets[r_idx + 1];
    m_data.insert(m_data.begin() + insert_data_pos, value);

    for (std::size_t i = r_idx + 1; i < m_row_offsets.size(); ++i) {
      m_row_offsets[i] += 1;
    }
  }

  void add_element(std::size_t r_idx, T&& value) {
    check_row_bounds(r_idx, "add_element");

    std::size_t insert_data_pos = m_row_offsets[r_idx + 1];
    m_data.insert(m_data.begin() + insert_data_pos, std::move(value));
    for (std::size_t i = r_idx + 1; i < m_row_offsets.size(); ++i) {
      m_row_offsets[i] += 1;
    }
  }

  void clear() noexcept {
    m_data.clear();
    m_row_offsets.assign(1, 0);
  }

  bool is_empty() const noexcept { return num_rows() == 0; }

  template <bool IsConstIter>
  class RowIterator {
   public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type =
        typename std::conditional<IsConstIter, std::span<const T>, std::span<T>>::type;
    using pointer = value_type*;
    using reference = value_type;

    using ParentArrayPtr =
        typename std::conditional<IsConstIter, const JaggedArray<T>*, JaggedArray<T>*>::type;

   private:
    ParentArrayPtr m_parent_array;
    std::size_t m_current_row_idx;

   public:
    RowIterator(ParentArrayPtr parent, std::size_t r_idx)
        : m_parent_array(parent), m_current_row_idx(r_idx) {}

    reference operator*() const {
      if (m_current_row_idx == m_parent_array->num_rows()) {
        if constexpr (IsConstIter) {
          return std::span<const T>();
        } else {
          return std::span<T>();
        }
      }
      // CORRECTED: Accessing members without trailing underscore
      std::size_t start_offset = m_parent_array->m_row_offsets[m_current_row_idx];
      std::size_t length = m_parent_array->m_row_offsets[m_current_row_idx + 1] - start_offset;

      if constexpr (IsConstIter) {
        return std::span<const T>(m_parent_array->m_data.data() + start_offset, length);
      } else {
        return std::span<T>(m_parent_array->m_data.data() + start_offset, length);
      }
    }

    RowIterator& operator++() {
      ++m_current_row_idx;
      return *this;
    }
    RowIterator operator++(int) {
      RowIterator tmp = *this;
      ++(*this);
      return tmp;
    }
    RowIterator& operator--() {
      --m_current_row_idx;
      return *this;
    }
    RowIterator operator--(int) {
      RowIterator tmp = *this;
      --(*this);
      return tmp;
    }

    RowIterator& operator+=(difference_type n) {
      m_current_row_idx += n;
      return *this;
    }
    RowIterator operator+(difference_type n) const {
      RowIterator temp = *this;
      temp += n;
      return temp;
    }
    friend RowIterator operator+(difference_type n, const RowIterator& iter) { return iter + n; }

    RowIterator& operator-=(difference_type n) {
      m_current_row_idx -= n;
      return *this;
    }
    RowIterator operator-(difference_type n) const {
      RowIterator temp = *this;
      temp -= n;
      return temp;
    }
    difference_type operator-(const RowIterator& other) const {
      return static_cast<difference_type>(m_current_row_idx) -
             static_cast<difference_type>(other.m_current_row_idx);
    }

    reference operator[](difference_type n) const { return *(*this + n); }

    bool operator==(const RowIterator& other) const {
      return m_parent_array == other.m_parent_array && m_current_row_idx == other.m_current_row_idx;
    }
    bool operator!=(const RowIterator& other) const { return !(*this == other); }
    bool operator<(const RowIterator& other) const {
      return m_current_row_idx < other.m_current_row_idx;
    }
    bool operator>(const RowIterator& other) const { return other < *this; }
    bool operator<=(const RowIterator& other) const { return !(other < *this); }
    bool operator>=(const RowIterator& other) const { return !(*this < other); }
  };

  using iterator = RowIterator<false>;
  using const_iterator = RowIterator<true>;

  iterator begin() { return iterator(this, 0); }
  iterator end() { return iterator(this, num_rows()); }

  const_iterator begin() const { return const_iterator(this, 0); }
  const_iterator end() const { return const_iterator(this, num_rows()); }
  const_iterator cbegin() const { return const_iterator(this, 0); }
  const_iterator cend() const { return const_iterator(this, num_rows()); }

  friend std::ostream& operator<< <>(std::ostream& os, const JaggedArray<T>& arr);
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const JaggedArray<T>& arr) {
  os << "JaggedArray (" << arr.num_rows() << " rows, " << arr.total_elements()
     << " total elements):\n";
  if (arr.is_empty()) {
    os << "  (empty)\n";
    return os;
  }
  for (std::size_t i = 0; i < arr.num_rows(); ++i) {
    os << "  [" << i << "] (" << arr.num_cols(i) << " elements): {";
    std::span<const T> current_row_span = arr[i];
    for (std::size_t j = 0; j < current_row_span.size(); ++j) {
      os << current_row_span[j];
      if (j < current_row_span.size() - 1) {
        os << ", ";
      }
    }
    os << "}\n";
  }
  return os;
}

}  // namespace oiseau::utils
