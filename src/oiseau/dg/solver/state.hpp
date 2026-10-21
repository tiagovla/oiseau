// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>

namespace oiseau::dg::solver {

/// Nodal values for an equation system, grouped by DGSpace element group.
///
/// Each array has shape [group cells, group nodes, number of fields]. The
/// component order matches `EquationSystem::fields()`. Arrays use contiguous
/// row-major storage so the solver can evaluate cell-local operators directly.
struct DGState {
  std::vector<xt::xarray<double>> groups;
};

/// Updates target = target_scale * target + source_scale * source.
/// Target and source may be the same state.
inline void scale_add_in_place(DGState& target, double target_scale, const DGState& source,
                               double source_scale) {
  for (std::size_t group = 0; group < target.groups.size(); ++group) {
    auto* values = target.groups[group].data();
    const auto* source_values = source.groups[group].data();
    for (std::size_t index = 0; index < target.groups[group].size(); ++index) {
      const auto old_value = values[index];
      const auto added_value = source_values[index];
      values[index] = target_scale * old_value + source_scale * added_value;
    }
  }
}

inline DGState add_scaled(const DGState& base, const DGState& increment, double scale) {
  auto result = base;
  for (std::size_t group = 0; group < result.groups.size(); ++group) {
    auto* values = result.groups[group].data();
    const auto* differences = increment.groups[group].data();
    for (std::size_t index = 0; index < result.groups[group].size(); ++index) {
      values[index] += scale * differences[index];
    }
  }
  return result;
}

}  // namespace oiseau::dg::solver
