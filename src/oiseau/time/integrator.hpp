// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <concepts>

namespace oiseau::time {

/// Interface for a time integrator acting on a system and its state.
///
/// The integrator updates state by one step starting at time. Concrete
/// integrators may keep scratch storage between calls. This compile-time
/// interface leaves the system and state types open to different solvers.
template <typename Integrator, typename System, typename State>
concept TimeIntegrator =
    requires(Integrator& integrator, System& system, State& state, double time, double time_step) {
      { integrator.step(system, state, time, time_step) } -> std::same_as<void>;
    };

}  // namespace oiseau::time
