// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>

#include "oiseau/dg/solver/rk4.hpp"
#include "oiseau/time/integrator.hpp"
#include "oiseau/time/lsrk54.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"

namespace {

struct ScalarState {
  double value = 0.0;
};

ScalarState add_scaled(const ScalarState& base, const ScalarState& increment, double scale) {
  return {base.value + scale * increment.value};
}

void scale_add_in_place(ScalarState& target, double target_scale, const ScalarState& source,
                        double source_scale) {
  const auto old_value = target.value;
  const auto added_value = source.value;
  target.value = target_scale * old_value + source_scale * added_value;
}

struct GrowthSystem {
  ScalarState rhs(const ScalarState& state) const { return {state.value}; }
};

struct TimeDependentSystem {
  ScalarState rhs(double time, const ScalarState&) const { return {time}; }
};

struct Euler {
  void step(const GrowthSystem& system, ScalarState& state, double, double time_step) {
    state = add_scaled(state, system.rhs(state), time_step);
  }
};

struct NoStep {};

struct RecordingSystem {
  std::vector<double> times;

  ScalarState rhs(double time, const ScalarState&) {
    times.push_back(time);
    return {1.0};
  }
};

struct DgDecaySystem {
  oiseau::dg::solver::DGState rhs(const oiseau::dg::solver::DGState& state) const {
    auto derivative = state;
    derivative.groups[0](0) = -state.groups[0](0);
    return derivative;
  }
};

static_assert(oiseau::time::TimeIntegrator<oiseau::time::RK4, GrowthSystem, ScalarState>);
static_assert(oiseau::time::TimeIntegrator<oiseau::time::LSRK54, GrowthSystem, ScalarState>);
static_assert(oiseau::time::TimeIntegrator<Euler, GrowthSystem, ScalarState>);
static_assert(!oiseau::time::TimeIntegrator<NoStep, GrowthSystem, ScalarState>);

}  // namespace

TEST(test_time_integrator, rk4_solves_autonomous_equation) {
  GrowthSystem system;
  ScalarState state{1.0};
  oiseau::time::RK4 integrator;

  const auto final_time = oiseau::time::run(integrator, system, state, 0.0, 0.1, 10);

  EXPECT_NEAR(final_time, 1.0, 1e-12);
  EXPECT_NEAR(state.value, std::exp(1.0), 3e-6);
}

TEST(test_time_integrator, rk4_passes_stage_times_to_rhs) {
  TimeDependentSystem system;
  ScalarState state;
  oiseau::time::RK4 integrator;

  integrator.step(system, state, 2.0, 0.5);

  EXPECT_DOUBLE_EQ(state.value, 1.125);
}

TEST(test_time_integrator, lsrk54_matches_nodaldg_stage_times) {
  RecordingSystem system;
  ScalarState state{3.0};
  oiseau::time::LSRK54 integrator;

  integrator.step(system, state, 2.0, 0.5);

  constexpr std::array<double, 5> stage_fractions = {
      0.0,
      1432997174477.0 / 9575080441755.0,
      2526269341429.0 / 6820363962896.0,
      2006345519317.0 / 3224310063776.0,
      2802321613138.0 / 2924317926251.0,
  };
  ASSERT_EQ(system.times.size(), stage_fractions.size());
  for (std::size_t stage = 0; stage < stage_fractions.size(); ++stage) {
    EXPECT_DOUBLE_EQ(system.times[stage], 2.0 + 0.5 * stage_fractions[stage]);
  }
  EXPECT_NEAR(state.value, 3.5, 1e-12);
}

TEST(test_time_integrator, lsrk54_integrates_time_dependent_rhs) {
  TimeDependentSystem system;
  ScalarState state;
  oiseau::time::LSRK54 integrator;

  integrator.step(system, state, 2.0, 0.5);

  EXPECT_NEAR(state.value, 1.125, 1e-12);
}

TEST(test_time_integrator, lsrk54_has_fourth_order_convergence) {
  const auto error = [](std::size_t step_count) {
    GrowthSystem system;
    ScalarState state{1.0};
    oiseau::time::LSRK54 integrator;
    (void)oiseau::time::run(integrator, system, state, 0.0, 1.0 / static_cast<double>(step_count),
                            step_count);
    return std::abs(state.value - std::exp(1.0));
  };

  const auto coarse_error = error(5);
  const auto fine_error = error(10);
  ASSERT_GT(coarse_error, 0.0);
  EXPECT_LT(fine_error, 0.08 * coarse_error);
}

TEST(test_time_integrator, shared_driver_accepts_another_scheme) {
  GrowthSystem system;
  ScalarState state{1.0};
  Euler integrator;
  std::ostringstream progress;

  const auto final_time = oiseau::time::run(integrator, system, state, 0.0, 0.1, 10, &progress);

  EXPECT_NEAR(final_time, 1.0, 1e-12);
  EXPECT_NEAR(state.value, std::pow(1.1, 10), 1e-12);
  EXPECT_NE(progress.str().find("100% (10/10)"), std::string::npos);
}

TEST(test_time_integrator, original_dg_rk4_api_remains_usable) {
  oiseau::dg::solver::DGState state;
  auto values = xt::xarray<double>::from_shape({1});
  values(0) = 1.0;
  state.groups.push_back(std::move(values));

  oiseau::dg::solver::RK4::run(DgDecaySystem{}, state, 0.1, 10);

  EXPECT_NEAR(state.groups[0](0), std::exp(-1.0), 3e-6);
}

TEST(test_time_integrator, lsrk54_integrates_dg_state) {
  DgDecaySystem system;
  oiseau::dg::solver::DGState state;
  auto values = xt::xarray<double>::from_shape({1});
  values(0) = 1.0;
  state.groups.push_back(std::move(values));
  oiseau::time::LSRK54 integrator;

  const auto final_time = oiseau::time::run(integrator, system, state, 0.0, 0.1, 10);

  EXPECT_NEAR(final_time, 1.0, 1e-12);
  EXPECT_NEAR(state.groups[0](0), std::exp(-1.0), 1e-5);
}
