// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

/// @file
/// @brief Integration rules for the reference interval and triangle.

namespace oiseau::utils::integration {

/**
 * @brief Return a cubature rule for the triangle with vertices (-1, -1), (1, -1), and (-1, 1).
 *
 * Orders 1 through 28 use tabulated rules. Higher orders use a mapped Gauss-Jacobi rule.
 *
 * @param order Requested total polynomial degree of exactness, greater than zero.
 * @return An N-by-2 array of (r, s) nodes and a length-N array of weights.
 * @throws std::out_of_range If @p order is zero or negative.
 */
std::pair<xt::xarray<double>, xt::xarray<double>> triangle_cubature(int order);

/**
 * @brief Return a Gauss-Legendre quadrature rule on [-1, 1].
 * @param order Requested polynomial degree of exactness, greater than zero.
 * @return One-dimensional arrays of nodes and their corresponding weights.
 * @throws std::out_of_range If @p order is zero or negative.
 */
std::pair<xt::xarray<double>, xt::xarray<double>> quadrature(int order);

/**
 * @brief Return Gauss-Jacobi nodes and weights on [-1, 1].
 * @param n Polynomial order. The rule contains @p n plus one node.
 * @param alpha Exponent of (1 - x) in the Jacobi weight, greater than -1.
 * @param beta Exponent of (1 + x) in the Jacobi weight, greater than -1.
 * @return One-dimensional arrays of nodes and their corresponding weights. For @p n equal to
 *         zero, the single weight is 2.
 */
std::pair<xt::xarray<double>, xt::xarray<double>> jacobi_gq(unsigned n, double alpha, double beta);

/**
 * @brief Return Gauss-Lobatto-Jacobi nodes on [-1, 1].
 * @param n Polynomial order, at least 1.
 * @param alpha Exponent of (1 - x) in the Jacobi weight, greater than -1.
 * @param beta Exponent of (1 + x) in the Jacobi weight, greater than -1.
 * @return A one-dimensional array of @p n plus one node, including both endpoints.
 */
xt::xarray<double> jacobi_gl(unsigned n, double alpha, double beta);

}  // namespace oiseau::utils::integration
