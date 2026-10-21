// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "oiseau/time/detail/rhs.hpp"

namespace oiseau::time {

/// Five-stage, fourth-order low-storage Runge-Kutta scheme used by nodaldg.
///
/// The state supplies scale_add_in_place(target, target_scale, source,
/// source_scale) in its namespace. The system supplies rhs(time, state) or
/// rhs(state). The in-place operation must support target and source aliasing.
/// Each step evaluates the RHS five times.
class LSRK54 {
 public:
  template <typename System, typename State>
  void step(System& system, State& state, double time, double time_step) {
    if (!std::isfinite(time)) {
      throw std::invalid_argument("Time must be finite");
    }
    if (!std::isfinite(time_step) || time_step <= 0.0) {
      throw std::invalid_argument("Time step must be finite and positive");
    }

    auto residual = detail::evaluate_rhs(system, time, state);
    scale_add_in_place(residual, 0.0, residual, time_step);
    scale_add_in_place(state, 1.0, residual, b[0]);

    for (std::size_t stage = 1; stage < a.size(); ++stage) {
      const auto derivative = detail::evaluate_rhs(system, time + c[stage] * time_step, state);
      scale_add_in_place(residual, a[stage], derivative, time_step);
      scale_add_in_place(state, 1.0, residual, b[stage]);
    }
  }

 private:
  // The a, b, and c coefficients match nodaldg.utils.integration.rk45.
  static constexpr std::array<double, 5> a = {
      0.0,
      -567301805773.0 / 1357537059087.0,
      -2404267990393.0 / 2016746695238.0,
      -3550918686646.0 / 2091501179385.0,
      -1275806237668.0 / 842570457699.0,
  };
  static constexpr std::array<double, 5> b = {
      1432997174477.0 / 9575080441755.0,  5161836677717.0 / 13612068292357.0,
      1720146321549.0 / 2090206949498.0,  3134564353537.0 / 4481467310338.0,
      2277821191437.0 / 14882151754819.0,
  };
  static constexpr std::array<double, 5> c = {
      0.0,
      1432997174477.0 / 9575080441755.0,
      2526269341429.0 / 6820363962896.0,
      2006345519317.0 / 3224310063776.0,
      2802321613138.0 / 2924317926251.0,
  };
};

}  // namespace oiseau::time
