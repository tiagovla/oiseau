// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/region.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace oiseau::mesh {

Region::Region(std::size_t tag, std::string name, std::vector<std::size_t> cells)
    : _tag(tag), _name(std::move(name)) {
  add_cells(std::move(cells));
}

std::size_t Region::tag() const { return _tag; }
const std::string& Region::name() const { return _name; }

void Region::add_cells(std::vector<std::size_t> cells) {
  if (_cells.empty()) {
    _cells = std::move(cells);
  } else {
    _cells.append_range(cells);
  }
  std::ranges::sort(_cells);
  _cells.erase(std::ranges::unique(_cells).begin(), _cells.end());
}

const std::vector<std::size_t>& Region::cells() const { return _cells; }

Region Region::operator|(const Region& other) const {
  if (_tag != other._tag || _name != other._name) {
    throw std::invalid_argument("Regions must have the same tag and name to merge");
  }

  Region merged(_tag, _name, {});
  std::ranges::set_union(_cells, other._cells, std::back_inserter(merged._cells));
  return merged;
}

}  // namespace oiseau::mesh
