// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace oiseau::mesh {

/** @brief Groups mesh cells under a numeric tag and a name. */
class Region {
 public:
  /**
   * @brief Constructs a region from mesh cell indices.
   * @param tag Numeric identifier for the region.
   * @param name Name of the region.
   * @param cells Initial cell indices. Duplicate indices are removed.
   */
  Region(std::size_t tag, std::string name, std::vector<std::size_t> cells);

  /** @brief Returns the numeric tag. */
  std::size_t tag() const;
  /** @brief Returns the region name. */
  const std::string& name() const;
  /**
   * @brief Adds cells to the region, keeping each index once.
   * @param cells Mesh cell indices to add, in any order.
   */
  void add_cells(std::vector<std::size_t> cells);
  /** @brief Returns the unique cell indices in ascending order. */
  const std::vector<std::size_t>& cells() const;

  /**
   * @brief Merges two regions with the same tag and name.
   * @param other Region to merge with this region.
   * @return A new region containing the union of both cell sets.
   * @throws std::invalid_argument If either the tag or the name differs.
   */
  Region operator|(const Region& other) const;

 private:
  std::size_t _tag;
  std::string _name;
  std::vector<std::size_t> _cells;
};

}  // namespace oiseau::mesh
