// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <chrono>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <stdexcept>

#include "oiseau/time/integrator.hpp"
#include "oiseau/time/progress_bar.hpp"

#ifndef NDEBUG
#include "oiseau/utils/logging.hpp"
#endif

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

#ifndef NDEBUG
  LOG_DEBUG("Starting time integration at t={} with dt={} for {} steps", initial_time, time_step,
            step_count);
  const auto start = std::chrono::steady_clock::now();
#endif

  ProgressBar bar(progress, step_count);
  bar.update(0);
  auto time = initial_time;
  for (std::size_t step = 0; step < step_count; ++step) {
    integrator.step(system, state, time, time_step);
    time += time_step;
    bar.update(step + 1);
  }

#ifndef NDEBUG
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  LOG_DEBUG("Finished time integration at t={} (wall time {:.3f} s)", time, elapsed.count());
#endif
  return time;
}

}  // namespace oiseau::time
