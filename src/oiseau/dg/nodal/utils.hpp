// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

namespace oiseau::dg::nodal::utils {

/**
 * @brief calculates the gauss-lobatto quadrature points and weights
 *
 * @param n order of the quadrature
 * @param alpha parameter of the jacobi polynomial
 * @param beta parameter of the jacobi polynomial
 */
std::pair<xt::xarray<double>, xt::xarray<double>> jacobi_gq(unsigned n, double alpha, double beta);

/**
 * @brief calculates the gauss-legendre quadrature points
 *
 * @param n order of the quadrature
 * @param alpha parameter of the jacobi polynomial
 * @param beta parameter of the jacobi polynomial
 */
xt::xarray<double> jacobi_gl(unsigned n, double alpha, double beta);

}  // namespace oiseau::dg::nodal::utils
