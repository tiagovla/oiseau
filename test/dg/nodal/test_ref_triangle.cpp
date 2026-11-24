// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_triangle.hpp"
#include "oiseau/test_macros.hpp"

TEST(test_ref_triangle, get_order) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefTriangle(order);
    EXPECT_EQ(ref.order(), order);
  }
}

TEST(test_ref_triangle, get_number_of_nodes) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefTriangle(order);
    unsigned n = (order + 1) * (order + 2) / 2;
    EXPECT_EQ(ref.number_of_nodes(), n);
  }
}

TEST(test_ref_triangle, order_two_nodes_are_equidistant) {
  const oiseau::dg::nodal::RefTriangle ref(2);
  const xt::xarray<double> expected = {{-1.0, -1.0}, {0.0, -1.0}, {1.0, -1.0},
                                       {-1.0, 0.0},  {0.0, 0.0},  {-1.0, 1.0}};
  EXPECT_FLOATS_NEARLY_EQ(ref.r(), expected, 1e-12);
}

TEST(test_ref_triangle, warp_factor) {
  xt::xarray<double> r = {0.1, 0.5, 0.9};
  xt::xarray<double> output = oiseau::dg::nodal::detail::warp_factor(3, r);
  xt::xarray<double> expected = {0.038, 0.192, 0.346};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_ref_triangle, basis_function) {
  xt::xarray<double> ab = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}};
  ab = xt::transpose(ab);
  xt::xarray<double> output = oiseau::dg::nodal::RefTriangle::basis_function(ab, 1, 2);
  xt::xarray<double> expected{0, -74.097064, -448.25662, -1502.7619, -3772.2142};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_ref_triangle, grad_basis_function) {
  const oiseau::dg::nodal::RefTriangle ref(3);
  const auto ab = oiseau::dg::nodal::detail::rs_to_ab(ref.r());
  const auto output = oiseau::dg::nodal::RefTriangle::grad_basis_function(ab, 1, 1);
  xt::xarray<double> expected = {
      {-2.1213, -2.1213, -2.1213, -2.1213, 0.8103, 1.4142, 0.8103, 5.5537, 5.5537, 8.4853},
      {-6.3640, -3.4324, 1.3110, 4.2426, -3.4324, 0.7071, 4.2426, 1.3110, 4.2426, 4.2426}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_ref_triangle, warped_equilateral_nodes) {
  unsigned order = 3;
  auto output = oiseau::dg::nodal::detail::generate_triangle_nodes(order);
  xt::xarray<double> expected = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -0.7236, -0.0000, 0.7236, -0.2764, 0.2764, 0},
      {-0.5774, -0.5774, -0.5774, -0.5774, -0.0986, -0.0000, -0.0986, 0.6760, 0.6760, 1.1547}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_ref_triangle, equilateral_xy_to_rs) {
  xt::xarray<double> xy = {{1, 2, 3}, {2, 3, 4}};
  xy = xt::transpose(xy);
  auto output = oiseau::dg::nodal::detail::equilateral_xy_to_rs(xy);
  xt::xarray<double> expected = {{-0.4880, 1.9761}, {-0.0654, 3.1308}, {0.3573, 4.2855}};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);
}

TEST(test_ref_triangle, rs_to_ab) {
  xt::xarray<double> rs = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -1.0000, -0.3333, 0.4472, -1.0000, -0.4472, -1.0000},
      {-1.0000, -1.0000, -1.0000, -1.0000, -0.4472, -0.3333, -0.4472, 0.4472, 0.4472, 1.0000}};
  rs = xt::transpose(rs);
  auto output = oiseau::dg::nodal::detail::rs_to_ab(rs);
  xt::xarray<double> expected = {
      {-1.0000, -0.4472, 0.4472, 1.0000, -1.0000, -0.0000, 1.0000, -1.0000, 1.0000, -1.0000},
      {-1.0000, -1.0000, -1.0000, -1.0000, -0.4472, -0.3333, -0.4472, 0.4472, 0.4472, 1.0000}};
  expected = xt::transpose(expected);
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.01);

  const double edge_distance = 1.0 / (1 << 24);
  const xt::xarray<double> near_apex = {{-1.0 + edge_distance / 4, 1.0 - edge_distance}};
  const auto collapsed = oiseau::dg::nodal::detail::rs_to_ab(near_apex);
  EXPECT_DOUBLE_EQ(collapsed(0, 0), -1.0);
  EXPECT_DOUBLE_EQ(collapsed(0, 1), near_apex(0, 1));
}

TEST(test_ref_triangle, vandermonde_interpolation) {
  for (unsigned order = 3; order <= 5; ++order) {
    // Coarse and fine reference elements
    auto coarse_triangle = oiseau::dg::nodal::RefTriangle(order);
    const oiseau::dg::nodal::RefElement& coarse_ref = coarse_triangle;
    auto fine_ref = oiseau::dg::nodal::RefTriangle(order + 2);

    // Vandermonde matrices at coarse and fine nodes
    auto V_coarse = coarse_ref.v();
    auto V_fine = coarse_ref.vandermonde(fine_ref.r());

    // Function at coarse nodes: f(r, s) = r² + s²
    auto r_coarse = xt::col(coarse_ref.r(), 0);
    auto s_coarse = xt::col(coarse_ref.r(), 1);
    auto f_coarse = xt::pow(r_coarse, 2) + xt::pow(s_coarse, 2);

    // Compute interpolation matrix: A = V_fine * V_coarse⁻¹ via solving transposed system
    auto AT = xt::linalg::solve(xt::transpose(V_coarse), xt::transpose(V_fine));
    auto interpolated = xt::linalg::dot(xt::transpose(AT), f_coarse);

    // Exact values at fine nodes
    auto r_fine = xt::col(fine_ref.r(), 0);
    auto s_fine = xt::col(fine_ref.r(), 1);
    auto f_fine_exact = xt::pow(r_fine, 2) + xt::pow(s_fine, 2);

    // Verify interpolation accuracy
    EXPECT_FLOATS_NEARLY_EQ(interpolated, f_fine_exact, 1e-10);
  }
}

TEST(test_ref_triangle, vandermonde_values) {
  const oiseau::dg::nodal::RefTriangle ref(3);
  const auto& output = ref.v();
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

TEST(test_ref_triangle, vandermonde_gradient_values) {
  const oiseau::dg::nodal::RefTriangle ref(3);
  const auto& output = ref.gv();
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

TEST(test_ref_triangle, differentiation_matrix_values) {
  const oiseau::dg::nodal::RefTriangle ref(3);
  const auto& output = ref.d();
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
