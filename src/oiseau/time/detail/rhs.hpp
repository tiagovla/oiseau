// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

namespace oiseau::time::detail {

template <typename System, typename State>
auto evaluate_rhs(System& system, double time, const State& state) {
  if constexpr (requires { system.rhs(time, state); }) {
    return system.rhs(time, state);
  } else {
    return system.rhs(state);
  }
}

}  // namespace oiseau::time::detail
