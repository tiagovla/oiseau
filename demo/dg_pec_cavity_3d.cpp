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
#include <vector>

#include "oiseau/utils/logging.hpp"

int main(int argc, char** argv) {
  oiseau::logging::init(argc, argv);
  using namespace oiseau::dg;
  using namespace oiseau::mesh;

  Mesh mesh{
      Topology({{0, 1, 2, 3, 4, 5, 6, 7}}, {CellKind::Hexahedron}),
      Geometry({
          {0, 0, 0},
          {1, 0, 0},
          {1, 1, 0},
          {0, 1, 0},
          {0, 0, 1},
          {1, 0, 1},
          {1, 1, 1},
          {0, 1, 1},
      }),
  };
  const DGSpace space(mesh, 4);

  const Field ex{"Ex"}, ey{"Ey"}, ez{"Ez"}, hx{"Hx"}, hy{"Hy"}, hz{"Hz"};

  // clang-format off
  const EquationSystem equations({ex, ey, ez, hx, hy, hz}, {
      dt(ex) == dy(hz) - dz(hy),
      dt(ey) == dz(hx) - dx(hz),
      dt(ez) == dx(hy) - dy(hx),
      dt(hx) == dz(ey) - dy(ez),
      dt(hy) == dx(ez) - dz(ex),
      dt(hz) == dy(ex) - dx(ey),
  });  // clang-format on
  const solver::LinearDG3D cavity(
      mesh, space, equations, [](const auto&) { return 1.0; },
      // Keep the normal electric component and reflect the tangential components.
      [](const auto& interior, const auto& normal) {
        const auto normal_e =
            interior[0] * normal[0] + interior[1] * normal[1] + interior[2] * normal[2];
        return std::vector<double>{2.0 * normal_e * normal[0] - interior[0],
                                   2.0 * normal_e * normal[1] - interior[1],
                                   2.0 * normal_e * normal[2] - interior[2],
                                   interior[3],
                                   interior[4],
                                   interior[5]};
      });

  const auto initial = cavity.initial_state(ez, [](double x, double y, double) {
    return std::sin(std::numbers::pi * x) * std::sin(std::numbers::pi * y);
  });

  const double half_period = 1.0 / std::sqrt(2.0);
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
