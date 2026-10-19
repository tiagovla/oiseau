// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

namespace oiseau::utils::integration {
std::pair<xt::xarray<double>, xt::xarray<double>> cubature(int order);
std::pair<xt::xarray<double>, xt::xarray<double>> quadrature(int order);
std::pair<xt::xarray<double>, xt::xarray<double>> jacobi_gq(unsigned n, double alpha, double beta);
}  // namespace oiseau::utils::integration
