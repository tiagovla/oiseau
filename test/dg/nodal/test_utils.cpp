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

TEST(test_dg_nodal_utils, d_matrix_3d_tensor_correctness) {
  const unsigned order = 4;

  // Generate tensor-product nodes (r, s, t)
  auto rst = generate_hexahedron_nodes(order);
  auto r = xt::col(rst, 0);
  auto s = xt::col(rst, 1);
  auto t = xt::col(rst, 2);

  // Build Vandermonde and gradient matrices
  auto V = vandermonde_3d_tensor(order, rst);
  auto gradV = grad_vandermonde_3d_tensor(order, rst);

  // Get differentiation matrices D_r, D_s, D_t
  auto D = d_matrix_3d(V, gradV);
  xt::xarray<double> Dr = xt::view(D, xt::all(), xt::all(), 0);
  xt::xarray<double> Ds = xt::view(D, xt::all(), xt::all(), 1);
  xt::xarray<double> Dt = xt::view(D, xt::all(), xt::all(), 2);

  // Test function: f(r, s, t) = r^3 + s^3 + t^3
  xt::xarray<double> f = xt::pow(r, 3) + xt::pow(s, 3) + xt::pow(t, 3);

  // Compute numerical derivatives
  auto df_dr = xt::linalg::dot(Dr, f);
  auto df_ds = xt::linalg::dot(Ds, f);
  auto df_dt = xt::linalg::dot(Dt, f);

  // Compute exact derivatives
  auto df_dr_expected = 3.0 * xt::pow(r, 2);
  auto df_ds_expected = 3.0 * xt::pow(s, 2);
  auto df_dt_expected = 3.0 * xt::pow(t, 2);

  // Check gradients with tight tolerance
  EXPECT_FLOATS_NEARLY_EQ(df_dr, df_dr_expected, 1e-10);
  EXPECT_FLOATS_NEARLY_EQ(df_ds, df_ds_expected, 1e-10);
  EXPECT_FLOATS_NEARLY_EQ(df_dt, df_dt_expected, 1e-10);
}
