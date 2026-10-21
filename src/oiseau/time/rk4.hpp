// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>

#include "oiseau/time/detail/rhs.hpp"

namespace oiseau::time {

/// One classical RK4 step with caller-provided RHS and state combination.
template <typename State, typename Rhs, typename AddScaled>
void rk4_step(State& state, double time, double time_step, Rhs&& rhs, AddScaled&& add_scaled_fn) {
  if (!std::isfinite(time)) {
    throw std::invalid_argument("Time must be finite");
  }
  if (!std::isfinite(time_step) || time_step <= 0.0) {
    throw std::invalid_argument("Time step must be finite and positive");
  }

  const auto k1 = rhs(time, state);
  const auto stage2 = add_scaled_fn(state, k1, 0.5 * time_step);
  const auto k2 = rhs(time + 0.5 * time_step, stage2);
  const auto stage3 = add_scaled_fn(state, k2, 0.5 * time_step);
  const auto k3 = rhs(time + 0.5 * time_step, stage3);
  const auto stage4 = add_scaled_fn(state, k3, time_step);
  const auto k4 = rhs(time + time_step, stage4);

  auto output = add_scaled_fn(state, k1, time_step / 6.0);
  output = add_scaled_fn(output, k2, time_step / 3.0);
  output = add_scaled_fn(output, k3, time_step / 3.0);
  output = add_scaled_fn(output, k4, time_step / 6.0);
  state = std::move(output);
}

/// Classical fourth-order Runge-Kutta integrator.
///
/// The system supplies rhs(time, state) or rhs(state). The state type supplies
/// add_scaled(base, increment, scale) in its namespace.
class RK4 {
 public:
  template <typename System, typename State>
  void step(System& system, State& state, double time, double time_step) {
    rk4_step(
        state, time, time_step,
        [&system](double stage_time, const State& value) {
          return detail::evaluate_rhs(system, stage_time, value);
        },
        [](const State& base, const State& increment, double scale) {
          return add_scaled(base, increment, scale);
        });
  }
};

}  // namespace oiseau::time
