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
#include <span>
#include <vector>

#include "oiseau/utils/logging.hpp"

std::vector<double> upwind(std::span<const double> interior, std::span<const double> exterior,
                           const oiseau::dg::solver::LinearDG<2>::Normal& normal,
                           std::span<const double> f_minus, std::span<const double> f_plus) {
  const auto [nx, ny] = normal;
  const double jump_ez = exterior[0] - interior[0];
  const double jump_ht = ny * (exterior[1] - interior[1]) - nx * (exterior[2] - interior[2]);
  return {0.5 * (f_minus[0] + f_plus[0] - jump_ez), 0.5 * (f_minus[1] + f_plus[1] - ny * jump_ht),
          0.5 * (f_minus[2] + f_plus[2] + nx * jump_ht)};
}

int main(int argc, char** argv) {
  oiseau::logging::init(argc, argv);
  using namespace oiseau::dg;
  using namespace oiseau::mesh;

  Mesh mesh{Topology({{0, 1, 4, 3}, {1, 2, 4}, {2, 5, 4}},
                     {CellKind::Quadrilateral, CellKind::Triangle, CellKind::Triangle}),
            oiseau::mesh::Geometry({{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}})};
  const Space space(mesh, 4);

  // TM Maxwell equations in vacuum with normalized permittivity and permeability.
  const Field ez{"Ez"}, hx{"Hx"}, hy{"Hy"};
  const EquationSystem equations(  // clang-format off
    {ez, hx, hy},
    {
        dt(ez) == dx(hy) - dy(hx),
        dt(hx) == -dy(ez),
        dt(hy) == dx(ez),
    }); //clang-format on
  const solver::LinearDG<2> cavity(
      mesh, space, equations, upwind,
      [](const auto&, const auto&, auto exterior) { exterior[0] *= -1.0; });

  const auto initial = cavity.initial_state(ez, [](double x, double y) {
    return std::sin(std::numbers::pi * x / 2.0) * std::sin(std::numbers::pi * y);
  });

  const double half_period = 2.0 / std::sqrt(5.0);
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
