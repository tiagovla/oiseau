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

TEST(test_dg_nodal_utils, test_vandermonde_1d) {
  unsigned n = 4;
  xt::xarray<double> r = {0.1, 0.3, 0.7, 0.9};
  auto output = vandermonde_1d(n, r);
  xt::xarray<double> expected = {{0.70710678, 0.12247449, -0.76685233, -0.27594723, 0.71687369},
                                 {0.70710678, 0.36742346, -0.57711567, -0.71559198, 0.1547238},
                                 {0.70710678, 0.85732141, 0.37156763, -0.36013452, -0.87411656},
                                 {0.70710678, 1.10227038, 1.13051426, 0.88396656, 0.44110205}};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}

TEST(test_dg_nodal_utils, test_grad_vandermonde_1d) {
  unsigned n = 4;
  xt::xarray<double> r = {0.1, 0.3, 0.7, 0.9};
  auto output = grad_vandermonde_1d(n, r);
  xt::xarray<double> expected = {{0., 1.22474487, 0.47434165, -2.66593089, -1.55386715},
                                 {0., 1.22474487, 1.42302495, -1.54343367, -3.77064691},
                                 {0., 1.22474487, 3.32039154, 4.06905241, 1.59629356},
                                 {0., 1.22474487, 4.26907484, 8.55904127, 12.74383196}};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}

TEST(test_dg_nodal_utils, test_generate_triangle_equidistant_nodes) {
  xt::xarray<double> output = generate_triangle_equidistant_nodes(2);
  xt::xarray<double> expected = {{-1., 0., 1., -0.5, 0.5, 0.},
                                 {-0.577, -0.577, -0.577, 0.289, 0.289, 1.155},
                                 {0., 0., 0., 0., 0., 0.}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_warp_factor) {
  xt::xarray<double> r = {0.1, 0.5, 0.9};
  xt::xarray<double> output = warp_factor(3, r);
  xt::xarray<double> expected = {0.038, 0.192, 0.346};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_simplexp_2d) {
  xt::xarray<double> ab = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
  ab = xt::transpose(ab);
  xt::xarray<double> output = simplexp_2d(ab, 1, 2);
  xt::xarray<double> expected{0, -74.097064, -448.25662, -1502.7619, -3772.2142};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_grad_simplexp_2d) {
  unsigned order = 3;
  auto rs = conversion_equilateral_xy_to_rs(generate_triangle_nodes(order));
  auto ab = conversion_rs_to_ab(rs);
  xt::xarray<double> output = grad_simplexp_2d(ab, 1, 1);
  xt::xarray<double> expected = {
      {-2.1213, -2.1213, -2.1213, -2.1213, 0.8103, 1.4142, 0.8103, 5.5537, 5.5537, 8.4853},
      {-6.3640, -3.4324, 1.3110, 4.2426, -3.4324, 0.7071, 4.2426, 1.3110, 4.2426, 4.2426}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_generate_triangle_nodes) {
  unsigned order = 3;
  auto output = generate_triangle_nodes(order);
  xt::xarray<double> expected = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -0.7236, -0.0000, 0.7236, -0.2764, 0.2764, 0},
      {-0.5774, -0.5774, -0.5774, -0.5774, -0.0986, -0.0000, -0.0986, 0.6760, 0.6760, 1.1547}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_conversion_equilateral_xy_to_rs) {
  xt::xarray<double> xy = {{1, 2, 3}, {2, 3, 4}};
  xy = xt::transpose(xy);
  auto output = conversion_equilateral_xy_to_rs(xy);
  xt::xarray<double> expected = {{-0.4880, 1.9761}, {-0.0654, 3.1308}, {0.3573, 4.2855}};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, d_matrix_2d_tensor_correctness) {
  const unsigned order = 5;

  // Generate nodal coordinates (r, s) and Vandermonde matrices
  auto rs = generate_quadrilateral_nodes(order);
  auto r = xt::col(rs, 0);
  auto s = xt::col(rs, 1);
  auto V = vandermonde_2d_tensor(order, rs);
  auto gradV = grad_vandermonde_2d_tensor(order, rs);

  // Compute differentiation matrices D_r and D_s
  auto D = d_matrix_2d(V, gradV);
  xt::xarray<double> Dr = xt::view(D, xt::all(), xt::all(), 0);
  xt::xarray<double> Ds = xt::view(D, xt::all(), xt::all(), 1);

  // Define test function: f(r,s) = r^3 + s^3
  xt::xarray<double> f = xt::pow(r, 3) + xt::pow(s, 3);

  // Compute numerical derivatives using D_r and D_s
  xt::xarray<double> df_dr = xt::linalg::dot(Dr, f);
  xt::xarray<double> df_ds = xt::linalg::dot(Ds, f);

  // Compute exact derivatives
  xt::xarray<double> df_dr_expected = 3.0 * xt::pow(r, 2);
  xt::xarray<double> df_ds_expected = 3.0 * xt::pow(s, 2);

  // Compare with tight tolerance
  EXPECT_FLOATS_NEARLY_EQ(df_dr, df_dr_expected, 1e-10);
  EXPECT_FLOATS_NEARLY_EQ(df_ds, df_ds_expected, 1e-10);
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

TEST(test_dg_nodal_utils, test_conversion_rs_to_ab) {
  xt::xarray<double> rs = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -1.0000, -0.3333, 0.4472, -1.0000, -0.4472, -1.0000},
      {-1.0000, -1.0000, -1.0000, -1.0000, -0.4472, -0.3333, -0.4472, 0.4472, 0.4472, 1.0000}};
  rs = xt::transpose(rs);
  auto output = conversion_rs_to_ab(rs);
  xt::xarray<double> expected = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -1.0000, -0.0000, 1.0000, -1.0000, 1.0000, -1.0000},
      {-1.0000, -1.0000, -1.0000, -1.0000, -0.4472, -0.3333, -0.4472, 0.4472, 0.4472, 1.0000}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_vandermonde_2d) {
  unsigned order = 3;
  auto rs = conversion_equilateral_xy_to_rs(generate_triangle_nodes(order));
  auto output = vandermonde_2d(order, rs);
  xt::xarray<double> expected = {
      {0.7071, -1.0000, 1.2247, -1.4142, -1.7321, 2.1213, -2.4495, 2.7386, -3.1623, -3.7417},
      {0.7071, -1.0000, 1.2247, -1.4142, -0.7746, 0.9487, -1.0954, -0.5477, 0.6325, 1.6733},
      {0.7071, -1.0000, 1.2247, -1.4142, 0.7746, -0.9487, 1.0954, -0.5477, 0.6325, -1.6733},
      {0.7071, -1.0000, 1.2247, -1.4142, 1.7321, -2.1213, 2.4495, 2.7386, -3.1623, 3.7417},
      {0.7071, -0.1708, -0.5477, 0.6325, -1.2533, -0.5863, 1.2628, 1.4340, 1.5478, -1.4177},
      {0.7071, -0.0000, -0.6804, 0.4190, -0.0000, -0.0000, 0.0000, -0.6086, -0.9370, 0.0000},
      {0.7071, -0.1708, -0.5477, 0.6325, 1.2533, 0.5863, -1.2628, 1.4340, 1.5478, 1.4177},
      {0.7071, 1.1708, 0.5477, -0.6325, -0.4787, -1.5350, -2.2426, 0.2092, 0.9821, -0.0790},
      {0.7071, 1.1708, 0.5477, -0.6325, 0.4787, 1.5350, 2.2426, 0.2092, 0.9821, 0.0790},
      {0.7071, 2.0000, 3.6742, 5.6569, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
  };
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_grad_vandermonde_2d) {
  unsigned order = 3;
  auto rs = conversion_equilateral_xy_to_rs(generate_triangle_nodes(order));
  auto output = grad_vandermonde_2d(order, rs);
  xt::xarray<double> expected_r = {
      {0, 0, 0, 0, 1.7321, -2.1213, 2.4495, -8.2158, 9.4868, 22.4499},
      {0, 0, 0, 0, 1.7321, -2.1213, 2.4495, -3.6742, 4.2426, 0.0000},
      {0, 0, 0, 0, 1.7321, -2.1213, 2.4495, 3.6742, -4.2426, -0.0000},
      {0, 0, 0, 0, 1.7321, -2.1213, 2.4495, 8.2158, -9.4868, 22.4499},
      {0, 0, 0, 0, 1.7321, 0.8103, -1.7452, -5.9450, -6.4168, 11.7549},
      {0, 0, 0, 0, 1.7321, 1.4142, -1.6330, -0.0000, -0.0000, -2.4944},
      {0, 0, 0, 0, 1.7321, 0.8103, -1.7452, 5.9450, 6.4168, 11.7549},
      {0, 0, 0, 0, 1.7321, 5.5537, 8.1138, -2.2708, -10.6595, 1.7150},
      {0, 0, 0, 0, 1.7321, 5.5537, 8.1138, 2.2708, 10.6595, 1.7150},
      {0, 0, 0, 0, 1.7321, 8.4853, 24.4949, 0, 0, 0},
  };
  xt::xarray<double> expected_s = {
      {0, 1.5000, -4.8990, 10.6066, 0.8660, -6.3640, 15.9217, -2.7386, 14.2302, 5.6125},
      {0, 1.5000, -4.8990, 10.6066, 0.8660, -3.4324, 7.7974, -0.4678, -1.6734, -2.5100},
      {0, 1.5000, -4.8990, 10.6066, 0.8660, 1.3110, -5.3479, 3.2064, -5.9161, 2.5100},
      {0, 1.5000, -4.8990, 10.6066, 0.8660, 4.2426, -13.4722, 5.4772, 4.7434, 16.8375},
      {0, 1.5000, -1.5139, -1.3110, 0.8660, -3.4324, -0.5256, -1.9817, 3.6563, 2.9387},
      {0, 1.5000, -0.8165, -2.3570, 0.8660, 0.7071, -0.8165, 0.9129, -1.0541, -1.2472},
      {0, 1.5000, -1.5139, -1.3110, 0.8660, 4.2426, -1.2195, 3.9634, 10.0732, 8.8162},
      {0, 1.5000, 3.9634, 3.4324, 0.8660, 1.3110, -2.1688, -0.7569, -2.7076, 0.4288},
      {0, 1.5000, 3.9634, 3.4324, 0.8660, 4.2426, 10.2826, 1.5139, 7.9518, 1.2863},
      {0, 1.5000, 7.3485, 21.2132, 0.8660, 4.2426, 12.2474, 0, 0, 0},
  };
  auto expected = xt::stack(xt::xtuple(expected_r, expected_s), 2);
  for (int i = 0; i < 2; ++i) {
    xt::xarray<double> c_output = xt::view(output, xt::all(), xt::all(), i);
    xt::xarray<double> c_expected = xt::view(expected, xt::all(), xt::all(), i);
    EXPECT_FLOATS_NEARLY_EQ(c_output, c_expected, 0.01);
  }
}

TEST(test_dg_nodal_utils, d_matrix_1d_consistency) {
  for (unsigned order = 2; order < 10; ++order) {
    auto r = jacobi_gl(order, 0, 0);
    auto vandermonde = vandermonde_1d(order, r);
    auto grad_vandermonde = grad_vandermonde_1d(order, r);
    auto d_matrix_from_v = d_matrix_1d(vandermonde, grad_vandermonde);
    auto d_matrix_direct = d_matrix_1d(order, r);
    EXPECT_FLOATS_NEARLY_EQ(d_matrix_from_v, d_matrix_direct, 1e-10);
  }
}

TEST(test_dg_nodal_utils, d_matrix_1d_differentiates_polynomials) {
  unsigned order = 10;
  auto r = jacobi_gl(order, 0, 0);
  auto d = d_matrix_1d(order, r);
  for (unsigned i = 1; i < order; ++i) {
    auto y = xt::pow(r, i);             // y = r^i
    auto dydr = xt::linalg::dot(d, y);  // dydr = i * r^(i-1)
    auto expected = (i)*xt::pow(r, (i - 1));
    EXPECT_FLOATS_NEARLY_EQ(dydr, expected, 1e-10);
  }
}

TEST(test_dg_nodal_utils, test_d_matrix_2d) {
  unsigned order = 3;
  auto rs = conversion_equilateral_xy_to_rs(generate_triangle_nodes(order));
  auto v = vandermonde_2d(order, rs);
  auto gv = grad_vandermonde_2d(order, rs);
  auto output = d_matrix_2d(v, gv);
  xt::xarray<double> expected_r = {
      {-3.0000, 4.0451, -1.5451, 0.5000, 0.0000, 0.0000, 0.0000, -0.0000, -0.0000, 0.0000},
      {-0.8090, 0.0000, 1.1180, -0.3090, 0.0000, -0.0000, -0.0000, 0.0000, 0.0000, -0.0000},
      {0.3090, -1.1180, 0.0000, 0.8090, 0.0000, -0.0000, 0.0000, 0.0000, -0.0000, -0.0000},
      {-0.5000, 1.5451, -4.0451, 3.0000, -0.0000, 0.0000, 0.0000, -0.0000, -0.0000, 0.0000},
      {-0.7090, 1.6180, -1.3090, 0.6000, -2.0000, 2.7000, -0.6180, -0.1910, -0.1910, 0.1000},
      {0.3333, -0.6211, 0.6211, -0.3333, -0.7272, -0.0000, 0.7272, -0.1061, 0.1061, -0.0000},
      {-0.6000, 1.3090, -1.6180, 0.7090, 0.6180, -2.7000, 2.0000, 0.1910, 0.1910, -0.1000},
      {0.4090, -0.1910, -0.6180, 0.6000, -1.3090, 2.7000, -1.3090, -2.0000, 1.6180, 0.1000},
      {-0.6000, 0.6180, 0.1910, -0.4090, 1.3090, -2.7000, 1.3090, -1.6180, 2.0000, -0.1000},
      {-0.5000, -0.0000, -0.0000, 0.5000, 1.5451, 0.0000, -1.5451, -4.0451, 4.0451, 0},
  };
  xt::xarray<double> expected_s = {
      {-3.0000, 0.0000, -0.0000, 0.0000, 4.0451, 0.0000, -0.0000, -1.5451, 0.0000, 0.5000},
      {-0.7090, -2.0000, -0.1910, 0.1000, 1.6180, 2.7000, -0.1910, -1.3090, -0.6180, 0.6000},
      {0.4090, -1.3090, -2.0000, 0.1000, -0.1910, 2.7000, 1.6180, -0.6180, -1.3090, 0.6000},
      {-0.5000, 1.5451, -4.0451, -0.0000, 0.0000, -0.0000, 4.0451, 0.0000, -1.5451, 0.5000},
      {-0.8090, 0.0000, 0.0000, 0.0000, 0.0000, -0.0000, 0.0000, 1.1180, -0.0000, -0.3090},
      {0.3333, -0.7272, -0.1061, -0.0000, -0.6211, -0.0000, 0.1061, 0.6211, 0.7272, -0.3333},
      {-0.6000, 1.3090, -1.6180, -0.1000, 0.6180, -2.7000, 2.0000, 0.1910, 1.3090, -0.4090},
      {0.3090, -0.0000, 0.0000, -0.0000, -1.1180, -0.0000, -0.0000, 0.0000, 0.0000, 0.8090},
      {-0.6000, 0.6180, 0.1910, -0.1000, 1.3090, -2.7000, 0.1910, -1.6180, 2.0000, 0.7090},
      {-0.5000, 0.0000, -0.0000, 0.0000, 1.5451, 0.0000, -0.0000, -4.0451, 0.0000, 3.0000},
  };
  auto expected = xt::stack(xt::xtuple(expected_r, expected_s), 2);
  for (int i = 0; i < 2; ++i) {
    xt::xarray<double> c_output = xt::view(output, xt::all(), xt::all(), i);
    xt::xarray<double> c_expected = xt::view(expected, xt::all(), xt::all(), i);
    EXPECT_FLOATS_NEARLY_EQ(c_output, c_expected, 0.01);
  }
}

TEST(test_dg_nodal_utils, test_generate_tetrahedron_equidistant_nodes) {
  unsigned order = 2;
  auto output = generate_tetrahedron_equidistant_nodes(order);
  xt::xarray<double> expected = {{-1, 0, 1, -1, 0, -1, -1, 0, -1, -1},
                                 {-1, -1, -1, 0, 0, 1, -1, -1, 0, -1},
                                 {-1, -1, -1, -1, -1, -1, 0, 0, 0, 1}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_eval_warp) {
  unsigned p = 3;
  xt::xarray<double> xnodes = {1.0000, 0.4472, -0.4472, -1.0000};
  xt::xarray<double> xout = {
      -1.0000, -0.3333, 0.3333, 1.0000,  -0.6667, -0.0000, 0.6667,  -0.3333, 0.3333,  0,
      -0.6667, -0.0000, 0.6667, -0.3333, 0.3333,  -0.0000, -0.3333, 0.3333,  -0.0000, 0,
  };
  xt::xarray<double> expected = {
      -3.8435e-01, -1.2812e-01, 1.2812e-01,  3.8435e-01,  -2.5623e-01, 6.9389e-17, 2.5623e-01,
      -1.2812e-01, 1.2812e-01,  8.3267e-17,  -2.5623e-01, 6.9389e-17,  2.5623e-01, -1.2812e-01,
      1.2812e-01,  6.9389e-17,  -1.2812e-01, 1.2812e-01,  6.9389e-17,  8.3267e-17,
  };
  auto output = eval_warp(p, xnodes, xout);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_eval_shift) {
  unsigned p = 3;
  unsigned pval = 0;
  xt::xarray<double> l = {
      {0, 1.0000, 0},           {0, 0.6667, 0.3333}, {0, 0.3333, 0.6667},
      {0, 0, 1.0000},           {0.3333, 0.6667, 0}, {0.3333, 0.3333, 0.3333},
      {0.3333, 0.0000, 0.6667}, {0.6667, 0.3333, 0}, {0.6667, 0.0000, 0.3333},
      {1.0000, 0, 0},           {0, 0.6667, 0},      {0, 0.3333, 0.3333},
      {0, 0.0000, 0.6667},      {0.3333, 0.3333, 0}, {0.3333, 0.0000, 0.3333},
      {0.6667, 0.0000, 0},      {0, 0.3333, 0},      {0, 0.0000, 0.3333},
      {0.3333, 0.0000, 0},      {0, 0, 0},
  };
  xt::xarray<double> expected = {
      {0, -0.1139, 0.1139, 0,       -0.0569, -0.0000, 0.0569, 0.0569, -0.0569, 0,
       0, 0.0000,  0.0000, -0.0000, -0.0000, 0.0000,  0,      0.0000, 0.0000,  0},
      {0, 0, 0, 0,         -0.098623, -0.000000, -0.098623, 0.098623, 0.098623, 0,
       0, 0, 0, -0.000000, 0.000000,  0.000000,  0,         0,        0.000000, 0}};
  expected = xt::transpose(expected);
  auto l1 = xt::col(l, 0);
  auto l2 = xt::col(l, 1);
  auto l3 = xt::col(l, 2);
  auto output = eval_shift(p, pval, l1, l2, l3);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_generate_tetrahedron_nodes) {
  xt::xarray<double> output_2 = generate_tetrahedron_nodes(2);
  xt::xarray<double> expected_2 = {
      {-1., -0.577, -0.408}, {0., -0.577, -0.408}, {1., -0.577, -0.408},  {-0.5, 0.289, -0.408},
      {0.5, 0.289, -0.408},  {0., 1.155, -0.408},  {-0.5, -0.289, 0.408}, {0.5, -0.289, 0.408},
      {0., 0.577, 0.408},    {0., 0., 1.225}};
  EXPECT_FLOATS_NEARLY_EQ(output_2, expected_2, 0.01);
  xt::xarray<double> output_3 = generate_tetrahedron_nodes(3);
  xt::xarray<double> expected_3 = {
      {-1.0000, -0.5774, -0.408248}, {-0.4472, -0.5774, -0.408248},
      {0.4472, -0.5774, -0.408248},  {1.0000, -0.5774, -0.408248},
      {-0.7236, -0.0986, -0.408248}, {-0.0000, -0.0000, -0.408248},
      {0.7236, -0.0986, -0.408248},  {-0.2764, 0.6760, -0.408248},
      {0.2764, 0.6760, -0.408248},   {0, 1.1547, -0.408248},
      {-0.7236, -0.4178, 0.043100},  {-0.0000, -0.3849, 0.136083},
      {0.7236, -0.4178, 0.043100},   {-0.3333, 0.1925, 0.136083},
      {0.3333, 0.1925, 0.136083},    {-0.0000, 0.8355, 0.043100},
      {-0.2764, -0.1596, 0.773397},  {0.2764, -0.1596, 0.773397},
      {-0.0000, 0.3192, 0.773397},   {0, 0, 1.224745},
  };
  EXPECT_FLOATS_NEARLY_EQ(output_3, expected_3, 0.01);
}

TEST(test_dg_nodal_utils, test_conversion_equilateral_xyz_to_rst) {
  xt::xarray<double> input = {
      {-1.0000, -0.5774, -0.408248}, {-0.4472, -0.5774, -0.408248},
      {0.4472, -0.5774, -0.408248},  {1.0000, -0.5774, -0.408248},
      {-0.7236, -0.0986, -0.408248}, {-0.0000, -0.0000, -0.408248},
      {0.7236, -0.0986, -0.408248},  {-0.2764, 0.6760, -0.408248},
      {0.2764, 0.6760, -0.408248},   {0, 1.1547, -0.408248},
      {-0.7236, -0.4178, 0.043100},  {-0.0000, -0.3849, 0.136083},
      {0.7236, -0.4178, 0.043100},   {-0.3333, 0.1925, 0.136083},
      {0.3333, 0.1925, 0.136083},    {-0.0000, 0.8355, 0.043100},
      {-0.2764, -0.1596, 0.773397},  {0.2764, -0.1596, 0.773397},
      {-0.0000, 0.3192, 0.773397},   {0, 0, 1.224745},
  };
  auto output = conversion_equilateral_xyz_to_rst(input);
  xt::xarray<double> expected = {
      {-1.0000, -1.0000, -1.0000}, {-0.4472, -1.0000, -1.0000}, {0.4472, -1.0000, -1.0000},
      {1.0000, -1.0000, -1.0000},  {-1.0000, -0.4472, -1.0000}, {-0.3333, -0.3333, -1.0000},
      {0.4472, -0.4472, -1.0000},  {-1.0000, 0.4472, -1.0000},  {-0.4472, 0.4472, -1.0000},
      {-1.0000, 1.0000, -1.0000},  {-1.0000, -1.0000, -0.4472}, {-0.3333, -1.0000, -0.3333},
      {0.4472, -1.0000, -0.4472},  {-1.0000, -0.3333, -0.3333}, {-0.3333, -0.3333, -0.3333},
      {-1.0000, 0.4472, -0.4472},  {-1.0000, -1.0000, 0.4472},  {-0.4472, -1.0000, 0.4472},
      {-1.0000, -0.4472, 0.4472},  {-1.0000, -1.0000, 1.0000},
  };
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.02);
}

TEST(test_dg_nodal_utils, test_simplexp_3d) {
  xt::xarray<double> input = {{1.0, 2.0, 2.0}};
  auto output = simplexp_3d(input, 1, 2, 3);
  xt::xarray<double> expected = {9575.4};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.05);
}

TEST(test_dg_nodal_utils, test_conversion_rst_to_abc) {
  xt::xarray<double> input = {
      {-1.0000, -1.0000, -1.0000}, {-0.4472, -1.0000, -1.0000}, {0.4472, -1.0000, -1.0000},
      {1.0000, -1.0000, -1.0000},  {-1.0000, -0.4472, -1.0000}, {-0.3333, -0.3333, -1.0000},
      {0.4472, -0.4472, -1.0000},  {-1.0000, 0.4472, -1.0000},  {-0.4472, 0.4472, -1.0000},
      {-1.0000, 1.0000, -1.0000},  {-1.0000, -1.0000, -0.4472}, {-0.3333, -1.0000, -0.3333},
      {0.4472, -1.0000, -0.4472},  {-1.0000, -0.3333, -0.3333}, {-0.3333, -0.3333, -0.3333},
      {-1.0000, 0.4472, -0.4472},  {-1.0000, -1.0000, 0.4472},  {-0.4472, -1.0000, 0.4472},
      {-1.0000, -0.4472, 0.4472},  {-1.0000, -1.0000, 1.0000},
  };
  auto output = conversion_rst_to_abc(input);
  xt::xarray<double> expected = {
      {-1.0000, -1.0000, -1.0000}, {-0.4472, -1.0000, -1.0000}, {0.4472, -1.0000, -1.0000},
      {1.0000, -1.0000, -1.0000},  {-1.0000, -0.4472, -1.0000}, {0.0000, -0.3333, -1.0000},
      {1.0000, -0.4472, -1.0000},  {-1.0000, 0.4472, -1.0000},  {1.0000, 0.4472, -1.0000},
      {-1.0000, 1.0000, -1.0000},  {-1.0000, -1.0000, -0.4472}, {0, -1.0000, -0.3333},
      {1.0000, -1.0000, -0.4472},  {-1.0000, 0, -0.3333},       {1.0000, 0, -0.3333},
      {-1.0000, 1.0000, -0.4472},  {-1.0000, -1.0000, 0.4472},  {1.0000, -1.0000, 0.4472},
      {-1.0000, 1.0000, 0.4472},   {-1.0000, -1.0000, 1.0000},
  };
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.02);
}

TEST(test_dg_nodal_utils, test_vandermonde_3d) {
  unsigned order = 3;
  auto rst = conversion_equilateral_xyz_to_rst(generate_tetrahedron_nodes(order));
  auto output = vandermonde_3d(order, rst);
  xt::xarray<double> expected = {
      {0.8660,  -1.1180, 1.3229,  -1.5000, -1.5811, 1.8708,  -2.1213, 2.2913,  -2.5981, -3.0000,
       -2.7386, 3.2404,  -3.6742, 3.9686,  -4.5000, -5.1962, 5.1235,  -5.8095, -6.7082, -7.9373},
      {0.8660,  -1.1180, 1.3229,  -1.5000, -1.5811, 1.8708,  -2.1213, 2.2913, -2.5981, -3.0000,
       -1.2247, 1.4491,  -1.6432, 1.7748,  -2.0125, -2.3238, -1.0247, 1.1619, 1.3416,  3.5496},
      {0.8660, -1.1180, 1.3229, -1.5000, -1.5811, 1.8708, -2.1213, 2.2913, -2.5981, -3.0000,
       1.2247, -1.4491, 1.6432, -1.7748, 2.0125,  2.3238, -1.0247, 1.1619, 1.3416,  -3.5496},
      {0.8660, -1.1180, 1.3229, -1.5000, -1.5811, 1.8708, -2.1213, 2.2913,  -2.5981, -3.0000,
       2.7386, -3.2404, 3.6742, -3.9686, 4.5000,  5.1962, 5.1235,  -5.8095, -6.7082, 7.9373},
      {0.8660,  -1.1180, 1.3229,  -1.5000, -0.2701, 0.3196, -0.3624, -1.0247, 1.1619, 1.3416,
       -1.9817, 2.3448,  -2.6587, -1.0969, 1.2438,  2.6788, 2.6827,  -3.0419, 3.2833, -3.0073},
      {0.8660, -1.1180, 1.3229, -1.5000, 0,       0,       0,       -1.2729, 1.4434,  0.8889,
       0.0000, -0.0000, 0.0000, 0.0000,  -0.0000, -0.0000, -1.1386, 1.2910,  -1.9876, -0.0000},
      {0.8660, -1.1180, 1.3229, -1.5000, -0.2701, 0.3196,  -0.3624, -1.0247, 1.1619, 1.3416,
       1.9817, -2.3448, 2.6587, 1.0969,  -1.2438, -2.6788, 2.6827,  -3.0419, 3.2833, 3.0073},
      {0.8660,  -1.1180, 1.3229,  -1.5000, 1.8512, -2.1904, 2.4837, 1.0247,  -1.1619, -1.3416,
       -0.7569, 0.8956,  -1.0155, -2.8717, 3.2562, -4.7573, 0.3914, -0.4438, 2.0833,  -0.1676},
      {0.8660, -1.1180, 1.3229, -1.5000, 1.8512,  -2.1904, 2.4837, 1.0247,  -1.1619, -1.3416,
       0.7569, -0.8956, 1.0155, 2.8717,  -3.2562, 4.7573,  0.3914, -0.4438, 2.0833,  0.1676},
      {0.8660, -1.1180, 1.3229, -1.5000, 3.1623,  -3.7417, 4.2426, 6.8739,  -7.7942, 12.0000,
       0.0000, -0.0000, 0.0000, 0.0000,  -0.0000, 0.0000,  0.0000, -0.0000, 0.0000,  0.0000},
      {0.8660,  0.1180,  -0.8176, 0.5171, -1.1441, -0.8913, 1.1213, 1.1997, 1.6476,  -1.1367,
       -1.9817, -1.5437, 1.9422,  2.0780, 2.8537,  -1.9687, 2.6827, 3.6842, -2.5416, -3.0073},
      {0.8660, 0.3727, -0.8819, 0.1111, -1.0541, -1.2472, 0.7857,  1.0184,  1.9245, -0.8889,
       0,      0,      0,       0,      0,       0,       -1.1386, -2.1517, 0.9938, 0},
      {0.8660, 0.1180, -0.8176, 0.5171,  -1.1441, -0.8913, 1.1213, 1.1997, 1.6476,  -1.1367,
       1.9817, 1.5437, -1.9422, -2.0780, -2.8537, 1.9687,  2.6827, 3.6842, -2.5416, 3.0073},
      {0.8660,  0.3727,  -0.8819, 0.1111,  0.5270,  0.6236,  -0.3928, -0.5092, -0.9623, -0.3333,
       -0.9129, -1.0801, 0.6804,  -1.3229, -2.5000, -0.1925, 0.5693,  1.0758,  1.2423,  -0.2940},
      {0.8660, 0.3727, -0.8819, 0.1111, 0.5270, 0.6236, -0.3928, -0.5092, -0.9623, -0.3333,
       0.9129, 1.0801, -0.6804, 1.3229, 2.5000, 0.1925, 0.5693,  1.0758,  1.2423,  0.2940},
      {0.8660, 0.1180, -0.8176, 0.5171, 2.2882, 1.7825, -2.2426, 3.5992, 4.9428, 4.5466,
       0,      0,      0,       0,      0,      0,      0,       0,      0,      0},
      {0.8660,  2.1180,  2.1405,  0.3829, -0.4370, -1.7279, -3.2426, 0.1750, 0.9505,  -0.0633,
       -0.7569, -2.9928, -5.6164, 0.3032, 1.6463,  -0.1097, 0.3914,  2.1253, -0.1416, -0.1676},
      {0.8660, 2.1180, 2.1405, 0.3829,  -0.4370, -1.7279, -3.2426, 0.1750, 0.9505,  -0.0633,
       0.7569, 2.9928, 5.6164, -0.3032, -1.6463, 0.1097,  0.3914,  2.1253, -0.1416, 0.1676},
      {0.8660, 2.1180, 2.1405, 0.3829, 0.8740, 3.4558, 6.4853, 0.5251, 2.8514, 0.2534,
       0,      0,      0,      0,      0,      0,      0,      0,      0,      0},
      {0.8660, 3.3541, 7.9373, 15.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  };
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.02);
}

TEST(test_dg_nodal_utils, test_grad_simplexp_3d) {
  unsigned order = 3;
  auto rst = conversion_equilateral_xyz_to_rst(generate_tetrahedron_nodes(order));
  auto abc = conversion_rst_to_abc(rst);
  xt::xarray<double> output = grad_simplexp_3d(abc, 1, 1, 3);
  xt::xarray<double> expected = {
      {5.408, 5.408, 5.408, 5.408, -2.066, -3.606,  -2.066, -14.159, -14.159, -21.633,
       2.614, 4.273, 2.614, -6.41, -6.41,  -10.455, -26.41, -26.41,  105.642, -0.},
      {16.225, 8.751, -3.343, -10.817, 8.751,  -1.803, -10.817, -3.343, -10.817, -10.817,
       7.841,  2.137, -5.227, 2.137,   -8.546, -5.227, -79.231, 52.821, 52.821,  -0.},
      {86.533,  40.194, -34.785, -81.125, -21.499, -1.803, 19.433, -65.035, 50.876, -10.817,
       -12.664, 2.137,  15.278,  3.873,   -10.282, -5.227, 1.848,  -28.258, 52.821, 0.}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_dg_nodal_utils, test_grad_vandermonde_3d) {
  unsigned order = 2;
  auto rst = conversion_equilateral_xyz_to_rst(generate_tetrahedron_nodes(order));
  auto output = grad_vandermonde_3d(order, rst);
  xt::xarray<double> expected_dr = {
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, -3.9686, -15.3704},
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, -3.9686, 0},
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, -3.9686, 15.3704},
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, 5.9529, -7.6852},
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, 5.9529, 7.6852},
      {0, 0, 0, 0, 0, 0, 2.7386, -3.2404, 15.8745, 0.0000},
      {0, 0, 0, 0, 0, 0, 2.7386, 6.4807, -1.9843, -7.6852},
      {0, 0, 0, 0, 0, 0, 2.7386, 6.4807, -1.9843, 7.6852},
      {0, 0, 0, 0, 0, 0, 2.7386, 6.4807, 7.9373, 0},
      {0, 0, 0, 0, 0, 0, 2.7386, 16.2019, 0, 0},
  };
  xt::xarray<double> expected_ds = {
      {0, 0, 0, 2.3717, -2.8062, -9.1652, 1.3693, -1.6202, -11.9059, -5.1235},
      {0, 0, 0, 2.3717, -2.8062, -9.1652, 1.3693, -1.6202, -1.9843, 2.5617},
      {0, 0, 0, 2.3717, -2.8062, -9.1652, 1.3693, -1.6202, 7.9373, 10.2470},
      {0, 0, 0, 2.3717, -2.8062, 2.2913, 1.3693, -1.6202, -1.9843, -2.5617},
      {0, 0, 0, 2.3717, -2.8062, 2.2913, 1.3693, -1.6202, 7.9373, 5.1235},
      {0, 0, 0, 2.3717, -2.8062, 13.7477, 1.3693, -1.6202, 7.9373, 0.0000},
      {0, 0, 0, 2.3717, 5.6125, -4.5826, 1.3693, 3.2404, -5.9529, -2.5617},
      {0, 0, 0, 2.3717, 5.6125, -4.5826, 1.3693, 3.2404, 3.9686, 5.1235},
      {0, 0, 0, 2.3717, 5.6125, 6.8739, 1.3693, 3.2404, 3.9686, 0},
      {0, 0, 0, 2.3717, 14.0312, 0, 1.3693, 8.1009, 0, 0},
  };
  xt::xarray<double> expected_dt = {
      {0, 2.2361, -6.6144, 0.7906, -6.5479, -2.2913, 1.3693, -11.3413, -3.9686, -5.1235},
      {0, 2.2361, -6.6144, 0.7906, -6.5479, -2.2913, 1.3693, -1.6202, -1.9843, 2.5617},
      {0, 2.2361, -6.6144, 0.7906, -6.5479, -2.2913, 1.3693, 8.1009, -0.0000, 10.2470},
      {0, 2.2361, -6.6144, 0.7906, 1.8708, 2.2913, 1.3693, -6.4807, 1.9843, -2.5617},
      {0, 2.2361, -6.6144, 0.7906, 1.8708, 2.2913, 1.3693, 3.2404, 3.9686, 5.1235},
      {0, 2.2361, -6.6144, 0.7906, 10.2896, 6.8739, 1.3693, -1.6202, 7.9373, 0.0000},
      {0, 2.2361, 3.3072, 0.7906, -0.9354, -1.1456, 1.3693, -1.6202, -1.9843, -2.5617},
      {0, 2.2361, 3.3072, 0.7906, -0.9354, -1.1456, 1.3693, 8.1009, -0.0000, 5.1235},
      {0, 2.2361, 3.3072, 0.7906, 7.4833, 3.4369, 1.3693, 3.2404, 3.9686, 0},
      {0, 2.2361, 13.2288, 0.7906, 4.6771, 0, 1.3693, 8.1009, 0, 0},
  };
  auto expected = xt::stack(xt::xtuple(expected_dr, expected_ds, expected_dt), 2);
  for (int i = 0; i < 3; ++i) {
    xt::xarray<double> c_output = xt::view(output, xt::all(), xt::all(), i);
    xt::xarray<double> c_expected = xt::view(expected, xt::all(), xt::all(), i);
    EXPECT_FLOATS_NEARLY_EQ(c_output, c_expected, 0.05);
  }
}

TEST(test_dg_nodal_utils, test_d_matrix_3d) {
  unsigned order = 2;
  auto rst = conversion_equilateral_xyz_to_rst(generate_tetrahedron_nodes(order));
  auto v = vandermonde_3d(order, rst);
  auto gv = grad_vandermonde_3d(order, rst);
  auto output = d_matrix_3d(v, gv);
  xt::xarray<double> expected_dr = {
      {-1.5000, 2.0000, -0.5000, -0.0000, 0.0000, -0.0000, -0.0000, 0.0000, 0, 0},
      {-0.5000, 0, 0.5000, 0, 0, 0, 0, 0, 0, 0},
      {0.5000, -2.0000, 1.5000, 0.0000, -0.0000, 0.0000, 0.0000, -0.0000, 0, 0},
      {-0.5000, 1.0000, -0.5000, -1.0000, 1.0000, -0.0000, 0.0000, -0.0000, 0, 0},
      {0.5000, -1.0000, 0.5000, -1.0000, 1.0000, -0.0000, 0.0000, -0.0000, 0, 0},
      {0.5000, 0.0000, -0.5000, -2.0000, 2.0000, -0.0000, 0.0000, -0.0000, 0, 0},
      {-0.5000, 1.0000, -0.5000, -0.0000, 0.0000, -0.0000, -1.0000, 1.0000, 0, 0},
      {0.5000, -1.0000, 0.5000, -0.0000, 0.0000, -0.0000, -1.0000, 1.0000, 0, 0},
      {0.5000, 0.0000, -0.5000, -1.0000, 1.0000, -0.0000, -1.0000, 1.0000, 0, 0},
      {0.5000, 0.0000, -0.5000, -0.0000, 0.0000, -0.0000, -2.0000, 2.0000, 0, 0},
  };
  xt::xarray<double> expected_ds = {
      {-1.5000, 0.0000, 0.0000, 2.0000, -0.0000, -0.5000, 0.0000, 0.0000, -0.0000, -0.0000},
      {-0.5000, -1.0000, 0.0000, 1.0000, 1.0000, -0.5000, 0.0000, -0.0000, -0.0000, -0.0000},
      {0.5000, -2.0000, -0.0000, 0.0000, 2.0000, -0.5000, 0.0000, -0.0000, -0.0000, 0.0000},
      {-0.5000, 0.0000, 0.0000, -0.0000, -0.0000, 0.5000, 0.0000, 0.0000, -0.0000, -0.0000},
      {0.5000, -1.0000, -0.0000, -1.0000, 1.0000, 0.5000, 0.0000, -0.0000, -0.0000, -0.0000},
      {0.5000, -0.0000, -0.0000, -2.0000, 0.0000, 1.5000, 0.0000, 0.0000, -0.0000, -0.0000},
      {-0.5000, 0.0000, 0.0000, 1.0000, -0.0000, -0.5000, -1.0000, -0.0000, 1.0000, -0.0000},
      {0.5000, -1.0000, -0.0000, 0.0000, 1.0000, -0.5000, -1.0000, 0.0000, 1.0000, -0.0000},
      {0.5000, 0.0000, 0.0000, -1.0000, 0.0000, 0.5000, -1.0000, -0.0000, 1.0000, -0.0000},
      {0.5000, 0.0000, -0.0000, -0.0000, 0.0000, -0.5000, -2.0000, -0.0000, 2.0000, 0.0000},
  };
  xt::xarray<double> expected_dt = {
      {-1.5000, -0.0000, -0.0000, -0.0000, -0.0000, 0.0000, 2.0000, 0.0000, -0.0000, -0.5000},
      {-0.5000, -1.0000, -0.0000, -0.0000, 0.0000, 0.0000, 1.0000, 1.0000, -0.0000, -0.5000},
      {0.5000, -2.0000, -0.0000, -0.0000, 0.0000, 0.0000, -0.0000, 2.0000, -0.0000, -0.5000},
      {-0.5000, -0.0000, -0.0000, -1.0000, 0.0000, -0.0000, 1.0000, 0.0000, 1.0000, -0.5000},
      {0.5000, -1.0000, -0.0000, -1.0000, 0.0000, -0.0000, -0.0000, 1.0000, 1.0000, -0.5000},
      {0.5000, -0.0000, -0.0000, -2.0000, 0.0000, -0.0000, 0, -0.0000, 2.0000, -0.5000},
      {-0.5000, 0, 0.0000, -0.0000, 0, 0.0000, 0.0000, 0.0000, 0.0000, 0.5000},
      {0.5000, -1.0000, -0.0000, -0.0000, 0.0000, 0.0000, -1.0000, 1.0000, 0.0000, 0.5000},
      {0.5000, 0.0000, -0.0000, -1.0000, 0.0000, -0.0000, -1.0000, -0.0000, 1.0000, 0.5000},
      {0.5000, 0.0000, -0.0000, 0.0000, 0.0000, 0.0000, -2.0000, 0.0000, 0, 1.5000},
  };
  auto expected = xt::stack(xt::xtuple(expected_dr, expected_ds, expected_dt), 2);
  for (int i = 0; i < 3; ++i) {
    xt::xarray<double> c_output = xt::view(output, xt::all(), xt::all(), i);
    xt::xarray<double> c_expected = xt::view(expected, xt::all(), xt::all(), i);
    EXPECT_FLOATS_NEARLY_EQ(c_output, c_expected, 0.05);
  }
}
