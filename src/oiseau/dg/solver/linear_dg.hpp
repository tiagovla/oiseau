// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <functional>
#include <memory>
#include <span>
#include <vector>

#include "oiseau/dg/dg_space.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/solver/state.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg::solver {

/// Explicit 2D DG solver for constant-coefficient linear first-order systems.
///
/// Equation declarations determine the volume fluxes. The numerical flux and
/// optional boundary ghost-state policy are supplied separately.
///
/// A `State` stores fields in the order returned by `equations().fields()`.
class LinearDG2D {
 public:
  using State = DGState;
  using Normal = std::array<double, 2>;
  using InitialCondition = std::function<double(double, double)>;
  using BoundaryState = std::function<std::vector<double>(std::span<const double>, const Normal&)>;
  using WaveSpeed = std::function<double(const Normal&)>;
  /// Arguments are left/right traces, outward normal, then each side's physical normal flux.
  using NumericalFlux = std::function<std::vector<double>(
      std::span<const double>, std::span<const double>, const Normal&, std::span<const double>,
      std::span<const double>)>;

  /// Uses the built-in Rusanov flux with this wave-speed bound.
  LinearDG2D(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem equations,
             WaveSpeed wave_speed, BoundaryState boundary_state = {});
  /// Uses the supplied numerical flux policy.
  LinearDG2D(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem equations,
             NumericalFlux numerical_flux, BoundaryState boundary_state = {});
  /// Existing-space variant using the built-in Rusanov flux.
  LinearDG2D(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations, WaveSpeed wave_speed,
             BoundaryState boundary_state = {});
  /// Existing-space variant using the supplied numerical flux policy.
  LinearDG2D(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
             NumericalFlux numerical_flux, BoundaryState boundary_state = {});
  ~LinearDG2D();

  LinearDG2D(const LinearDG2D&) = delete;
  LinearDG2D& operator=(const LinearDG2D&) = delete;
  LinearDG2D(LinearDG2D&&) noexcept;
  LinearDG2D& operator=(LinearDG2D&&) noexcept;

  const DGSpace& space() const;
  const EquationSystem& equations() const;
  /// Builds the local Lax-Friedrichs (Rusanov) numerical flux policy.
  static NumericalFlux rusanov(WaveSpeed wave_speed);
  State zero_state() const;
  /// Creates a zero state and samples this scalar field at every DG node.
  /// Other fields remain zero.
  State initial_state(const Field& field, const InitialCondition& value) const;
  /// Sets one scalar field by sampling its initial condition at every DG node.
  void set_initial_condition(State& state, const Field& field, const InitialCondition& value) const;
  State rhs(const State& state) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace oiseau::dg::solver
