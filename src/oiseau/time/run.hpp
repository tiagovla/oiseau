// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cmath>
#include <cstddef>
#include <ostream>
#include <stdexcept>

#include "oiseau/time/integrator.hpp"
#include "oiseau/time/progress_bar.hpp"

namespace oiseau::time {

/// Runs a fixed number of equal time steps and returns the final time.
template <typename Integrator, typename System, typename State>
  requires TimeIntegrator<Integrator, System, State>
double run(Integrator& integrator, System& system, State& state, double initial_time,
           double time_step, std::size_t step_count, std::ostream* progress = nullptr) {
  if (!std::isfinite(initial_time)) {
    throw std::invalid_argument("Initial time must be finite");
  }
  if (!std::isfinite(time_step) || time_step <= 0.0) {
    throw std::invalid_argument("Time step must be finite and positive");
  }

  ProgressBar bar(progress, step_count);
  bar.update(0);
  auto time = initial_time;
  for (std::size_t step = 0; step < step_count; ++step) {
    integrator.step(system, state, time, time_step);
    time += time_step;
    bar.update(step + 1);
  }
  return time;
}

}  // namespace oiseau::time
