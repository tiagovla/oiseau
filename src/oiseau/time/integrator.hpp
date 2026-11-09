// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <concepts>

namespace oiseau::time {

template <typename State>
concept AddScaledState = requires(const State& base, const State& increment, double scale) {
  { add_scaled(base, increment, scale) } -> std::same_as<State>;
};

template <typename State>
concept ScaleAddInPlaceState =
    requires(State& target, const State& source, double target_scale, double source_scale) {
      { scale_add_in_place(target, target_scale, source, source_scale) } -> std::same_as<void>;
    };

template <typename Integrator, typename System, typename State>
concept TimeIntegrator =
    requires(Integrator& integrator, System& system, State& state, double time, double time_step) {
      { integrator.step(system, state, time, time_step) } -> std::same_as<void>;
    };

}  // namespace oiseau::time
