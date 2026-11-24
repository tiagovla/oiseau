// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

#include "oiseau/dg/dg_space.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/solver/state.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg::solver {

/// Explicit strong-form DG solver for constant-coefficient linear first-order systems.
///
/// Equation declarations determine the volume fluxes. The numerical flux and
/// optional boundary ghost-state policy are supplied separately.
///
/// A `State` stores fields in the order returned by `equations().fields()`.
template <std::size_t Dimension>
class LinearDG {
  static_assert(Dimension >= 1 && Dimension <= 3);

 public:
  using State = DGState;
  using Normal = std::array<double, Dimension>;
  using InitialCondition =
      std::conditional_t<Dimension == 1, std::function<double(double)>,
                         std::conditional_t<Dimension == 2, std::function<double(double, double)>,
                                            std::function<double(double, double, double)>>>;
  using BoundaryState = std::function<std::vector<double>(std::span<const double>, const Normal&)>;
  using WaveSpeed = std::function<double(const Normal&)>;
  /// Arguments are left/right traces, outward normal, then each side's physical normal flux.
  using NumericalFlux = std::function<std::vector<double>(
      std::span<const double>, std::span<const double>, const Normal&, std::span<const double>,
      std::span<const double>)>;

  /// Uses the built-in Rusanov flux with this wave-speed bound.
  LinearDG(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem equations,
           WaveSpeed wave_speed, BoundaryState boundary_state = {});
  /// Uses the supplied numerical flux policy.
  LinearDG(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem equations,
           NumericalFlux numerical_flux, BoundaryState boundary_state = {});
  /// Existing-space variant using the built-in Rusanov flux.
  LinearDG(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations, WaveSpeed wave_speed,
           BoundaryState boundary_state = {});
  /// Existing-space variant using the supplied numerical flux policy.
  LinearDG(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
           NumericalFlux numerical_flux, BoundaryState boundary_state = {});
  ~LinearDG();

  LinearDG(const LinearDG&) = delete;
  LinearDG& operator=(const LinearDG&) = delete;
  LinearDG(LinearDG&&) noexcept;
  LinearDG& operator=(LinearDG&&) noexcept;

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

/// Solver for intervals.
using LinearDG1D = LinearDG<1>;
/// Solver for affine triangles and parallelograms.
using LinearDG2D = LinearDG<2>;
/// Solver for affine tetrahedra and parallelepipeds.
using LinearDG3D = LinearDG<3>;

}  // namespace oiseau::dg::solver
