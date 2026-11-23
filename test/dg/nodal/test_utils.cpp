// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xstrided_view.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/utils.hpp"
#include "oiseau/test_macros.hpp"
#include "xtensor/core/xmath.hpp"

using namespace oiseau::dg::nodal::utils;

TEST(test_dg_nodal_utils, test_jacobi_gq) {
  unsigned n = 3;
  double alpha = 1.0, beta = 2.0;
  auto [x, w] = jacobi_gq(n, alpha, beta);
  xt::xarray<double> x_expected = {-0.65077886, -0.15637043, 0.37348938, 0.79729627};
  xt::xarray<double> w_expected = {0.08666291, 0.44123335, 0.59015336, 0.21528371};
  EXPECT_FLOATS_NEARLY_EQ(x, x_expected, 0.0001);
  EXPECT_FLOATS_NEARLY_EQ(w, w_expected, 0.0001);
}

TEST(test_dg_nodal_utils, test_jacobi_gl_n1) {
  unsigned n = 1;
  double alpha = 0.0, beta = 0.0;
  auto output = jacobi_gl(n, alpha, beta);
  xt::xarray<double> expected = {-1., 1.};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 1e-10);
}

TEST(test_dg_nodal_utils, test_jacobi_gl) {
  unsigned n = 4;
  double alpha = 1.0, beta = 2.0;
  auto output = jacobi_gl(n, alpha, beta);
  xt::xarray<double> expected = {-1., -0.4611587, 0.10438533, 0.62950065, 1.};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}
