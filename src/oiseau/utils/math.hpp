// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <ranges>
#include <xtensor/containers/xarray.hpp>

namespace oiseau::utils {

/// Concept that checks for contiguous floating-point containers.
template <class Container>
concept FloatingArrayLike = std::ranges::contiguous_range<Container> &&
                            std::is_floating_point_v<std::ranges::range_value_t<Container>>;

/// Computes the Jacobi polynomial P_n^{(alpha,beta)}(x) at a single point.
template <std::floating_point Real>
Real jacobi(unsigned n, Real alpha, Real beta, Real x) {
  if (n == 0) {
    return Real(1);
  }

  Real y0 = 1;
  Real y1 = (alpha + 1) + (alpha + beta + 2) * (x - 1) / Real(2);

  Real yk = y1;
  unsigned k = 2;
  Real k_max = n * (1 + std::numeric_limits<Real>::epsilon());

  while (k < k_max) {
    Real denom = 2 * k * (k + alpha + beta) * (2 * k + alpha + beta - 2);
    Real gamma1 =
        (2 * k + alpha + beta - 1) *
        ((2 * k + alpha + beta) * (2 * k + alpha + beta - 2) * x + alpha * alpha - beta * beta);
    Real gamma0 = -2 * (k + alpha - 1) * (k + beta - 1) * (2 * k + alpha + beta);
    yk = (gamma1 * y1 + gamma0 * y0) / denom;
    y0 = y1;
    y1 = yk;
    ++k;
  }
  return yk;
}

/// Computes the normalized Jacobi polynomial at each entry of a container.
template <std::floating_point Real, FloatingArrayLike Container>
Container jacobi_p(unsigned n, Real alpha, Real beta, const Container& input) {
  Container v = input;
  Real norm = std::pow(2, alpha + beta + 1) / (2 * n + alpha + beta + 1);
  norm *= std::tgamma(n + alpha + 1) * std::tgamma(n + beta + 1);
  norm /= (std::tgamma(n + 1) * std::tgamma(n + alpha + beta + 1));
  auto partial = [n, alpha, beta, norm](auto& k) {
    k = jacobi(n, alpha, beta, k) / std::sqrt(norm);
  };
  std::for_each(v.begin(), v.end(), partial);
  return v;
}

/// Computes the gradient of the normalized Jacobi polynomial.
template <std::floating_point Real, FloatingArrayLike Container>
Container grad_jacobi_p(int n, Real alpha, Real beta, const Container& input) {
  Container v = input;
  if (n == 0) {
    std::fill(v.begin(), v.end(), Real(0));
    return v;
  }
  Real norm = std::sqrt(n * (n + alpha + beta + 1));
  Container output = jacobi_p(n - 1, alpha + 1, beta + 1, v);
  std::for_each(output.begin(), output.end(), [norm](auto& k) { k *= norm; });
  return output;
}

}  // namespace oiseau::utils
