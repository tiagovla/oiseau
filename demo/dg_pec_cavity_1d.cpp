// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <oiseau.hpp>

#include "oiseau/utils/logging.hpp"

int main(int argc, char** argv) {
  oiseau::logging::init(argc, argv);
  using namespace oiseau::dg;
  using namespace oiseau::mesh;

  Mesh mesh{Topology({{0, 1}, {1, 2}}, {CellKind::Interval, CellKind::Interval}),
            Geometry({{0}, {1}, {2}})};
  const DGSpace space(mesh, {3, 4});

  const Field ez{"Ez"}, hy{"Hy"};
  const EquationSystem equations({ez, hy}, {dt(ez) == dx(hy), dt(hy) == dx(ez)});
  const solver::LinearDG<1> cavity(
      mesh, space, equations, [](const auto&) { return 1.0; },
      // Reflect the tangential electric field at both PEC endpoints.
      [](const auto&, const auto&, auto exterior) { exterior[0] *= -1.0; });

  const auto initial =
      cavity.initial_state(ez, [](double x) { return std::sin(std::numbers::pi * x / 2.0); });

  constexpr double half_period = 2.0;
  constexpr std::size_t steps = 256;
  oiseau::time::RK4 rk4;
  auto state = initial;
  oiseau::time::run(rk4, cavity, state, 0.0, half_period / steps, steps);

  // The exact electric field changes sign after half a period.
  const auto field = equations.field_index(ez);
  double max_error = 0.0;
  for (std::size_t group = 0; group < state.groups.size(); ++group) {
    const auto& values = state.groups[group];
    const auto& start = initial.groups[group];
    for (std::size_t cell = 0; cell < values.shape()[0]; ++cell) {
      for (std::size_t node = 0; node < values.shape()[1]; ++node) {
        const auto error = std::abs(values(cell, node, field) + start(cell, node, field));
        if (!std::isfinite(error)) {
          std::cerr << "Nonfinite electric field\n";
          return 1;
        }
        max_error = std::max(max_error, error);
      }
    }
  }
  std::cout << "Maximum Ez error after half a period: " << max_error << '\n';
  return max_error < 0.03 ? 0 : 1;
}
