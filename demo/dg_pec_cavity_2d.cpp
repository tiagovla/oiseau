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
#include <vector>

#include "oiseau/dg/dg_space.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"

int main() {
  using namespace oiseau::dg;
  using namespace oiseau::mesh;

  const auto quad = get_cell_type(CellKind::Quadrilateral);
  const auto triangle = get_cell_type(CellKind::Triangle);
  Mesh mesh{Topology({{0, 1, 4, 3}, {1, 2, 4}, {2, 5, 4}}, {quad, triangle, triangle}),
            Geometry({0, 0, 1, 0, 2, 0, 0, 1, 1, 1, 2, 1}, 2)};
  const DGSpace space(mesh, 4);

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
      mesh, space, equations, [](const auto&) { return 1.0; },
      // Ez is tangential to every wall. Its odd reflection enforces the PEC condition.
      [](const auto& interior, const auto&) {
        return std::vector<double>{-interior[0], interior[1], interior[2]};
      });

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
