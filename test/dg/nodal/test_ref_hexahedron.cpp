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
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_hexahedron.hpp"
#include "oiseau/test_macros.hpp"

namespace oiseau {
namespace {

TEST(test_ref_hexahedron, get_order) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = dg::nodal::RefHexahedron(order);
    EXPECT_EQ(ref.order(), order);
  }
}

TEST(test_ref_hexahedron, get_number_of_nodes) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = dg::nodal::RefHexahedron(order);
    unsigned n = (order + 1) * (order + 1) * (order + 1);
    EXPECT_EQ(ref.number_of_nodes(), n);
  }
}

TEST(test_ref_hexahedron, vandermonde_interpolation) {
  for (unsigned order = 2; order <= 3; ++order) {
    // Coarse and fine reference elements
    const dg::nodal::RefHexahedron coarse_element(order);
    const dg::nodal::RefElement& coarse_ref = coarse_element;
    auto fine_ref = dg::nodal::RefHexahedron(order + 2);

    // Vandermonde matrices at coarse and fine nodes
    auto V_coarse = coarse_ref.v();
    auto V_fine = coarse_ref.vandermonde(fine_ref.r());

    // Function at coarse nodes: f(r, s) = r² + s² + t²
    auto r_coarse = xt::col(coarse_ref.r(), 0);
    auto s_coarse = xt::col(coarse_ref.r(), 1);
    auto t_coarse = xt::col(coarse_ref.r(), 2);
    auto f_coarse = xt::pow(r_coarse, 2) + xt::pow(s_coarse, 2) + xt::pow(t_coarse, 2);

    // Compute interpolation matrix: A = V_fine * V_coarse⁻¹ via solving transposed system
    auto AT = xt::linalg::solve(xt::transpose(V_coarse), xt::transpose(V_fine));
    auto interpolated = xt::linalg::dot(xt::transpose(AT), f_coarse);

    // Exact values at fine nodes
    auto r_fine = xt::col(fine_ref.r(), 0);
    auto s_fine = xt::col(fine_ref.r(), 1);
    auto t_fine = xt::col(fine_ref.r(), 2);
    auto f_fine_exact = xt::pow(r_fine, 2) + xt::pow(s_fine, 2) + xt::pow(t_fine, 2);

    // Verify interpolation accuracy
    EXPECT_FLOATS_NEARLY_EQ(interpolated, f_fine_exact, 1e-10);
  }
}

TEST(test_ref_hexahedron, differentiation_matrix_differentiates_polynomials) {
  const dg::nodal::RefHexahedron ref(4);
  auto r = xt::col(ref.r(), 0);
  auto s = xt::col(ref.r(), 1);
  auto t = xt::col(ref.r(), 2);

  const auto& D = ref.d();
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

}  // namespace
}  // namespace oiseau
