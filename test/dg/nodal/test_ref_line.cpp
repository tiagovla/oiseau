// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_line.hpp"
#include "oiseau/test_macros.hpp"

TEST(test_ref_line, get_order) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefLine(order);
    EXPECT_EQ(ref.order(), order);
  }
}

TEST(test_ref_line, get_number_of_nodes) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefLine(order);
    EXPECT_EQ(ref.number_of_nodes(), order + 1);
  }
}

TEST(test_ref_line, differentiation_matrix_differentiates_polynomials) {
  for (unsigned order = 1; order <= 10; ++order) {
    SCOPED_TRACE(order);
    const oiseau::dg::nodal::RefLine ref(order);
    const auto& derivative = ref.d();
    ASSERT_EQ(derivative.dimension(), 2);
    ASSERT_EQ(derivative.shape()[0], order + 1);
    ASSERT_EQ(derivative.shape()[1], order + 1);

    for (unsigned degree = 0; degree <= order; ++degree) {
      SCOPED_TRACE(degree);
      const xt::xarray<double> values = xt::pow(ref.r(), degree);
      const xt::xarray<double> differentiated = xt::linalg::dot(derivative, values);
      for (std::size_t node = 0; node < ref.r().size(); ++node) {
        const auto expected = degree == 0 ? 0.0 : degree * std::pow(ref.r()(node), degree - 1);
        EXPECT_NEAR(differentiated(node), expected, 1e-11);
      }
    }
  }
}

TEST(test_ref_line, vandermonde_interpolation) {
  for (unsigned order = 3; order <= 5; ++order) {
    // Coarse and fine reference elements
    const oiseau::dg::nodal::RefLine coarse_element(order);
    const oiseau::dg::nodal::RefElement& coarse_ref = coarse_element;
    auto fine_ref = oiseau::dg::nodal::RefLine(order + 2);

    // Vandermonde matrices at coarse and fine nodes
    auto V_coarse = coarse_ref.v();
    auto V_fine = coarse_ref.vandermonde(fine_ref.r());

    // Function at coarse nodes: f(r, s) = r²
    auto r_coarse = coarse_ref.r();
    auto f_coarse = xt::pow(r_coarse, 2);

    // Compute interpolation matrix: A = V_fine * V_coarse⁻¹ via solving transposed system
    auto AT = xt::linalg::solve(xt::transpose(V_coarse), xt::transpose(V_fine));
    auto interpolated = xt::linalg::dot(xt::transpose(AT), f_coarse);

    // Exact values at fine nodes
    auto r_fine = fine_ref.r();
    auto f_fine_exact = xt::pow(r_fine, 2);

    // Verify interpolation accuracy
    EXPECT_FLOATS_NEARLY_EQ(interpolated, f_fine_exact, 1e-10);
  }
}

TEST(test_ref_line, vandermonde_values) {
  unsigned n = 4;
  xt::xarray<double> r = {0.1, 0.3, 0.7, 0.9};
  const oiseau::dg::nodal::RefLine ref(n);
  const oiseau::dg::nodal::RefElement& element = ref;
  const auto output = element.vandermonde(r);
  xt::xarray<double> expected = {{0.70710678, 0.12247449, -0.76685233, -0.27594723, 0.71687369},
                                 {0.70710678, 0.36742346, -0.57711567, -0.71559198, 0.1547238},
                                 {0.70710678, 0.85732141, 0.37156763, -0.36013452, -0.87411656},
                                 {0.70710678, 1.10227038, 1.13051426, 0.88396656, 0.44110205}};
  EXPECT_FLOATS_NEARLY_EQ(output, expected, 0.0001);
}

TEST(test_ref_line, basis_gradient_values) {
  unsigned n = 4;
  xt::xarray<double> r = {0.1, 0.3, 0.7, 0.9};
  xt::xarray<double> expected = {{0., 1.22474487, 0.47434165, -2.66593089, -1.55386715},
                                 {0., 1.22474487, 1.42302495, -1.54343367, -3.77064691},
                                 {0., 1.22474487, 3.32039154, 4.06905241, 1.59629356},
                                 {0., 1.22474487, 4.26907484, 8.55904127, 12.74383196}};
  for (unsigned degree = 0; degree <= n; ++degree) {
    SCOPED_TRACE(degree);
    const auto output = oiseau::dg::nodal::RefLine::grad_basis_function(r, degree);
    const xt::xarray<double> expected_column = xt::col(expected, degree);
    EXPECT_FLOATS_NEARLY_EQ(output, expected_column, 0.0001);
  }
}

TEST(test_ref_line, differentiation_matrix_reproduces_basis_gradients) {
  for (unsigned order = 2; order < 10; ++order) {
    SCOPED_TRACE(order);
    const oiseau::dg::nodal::RefLine ref(order);
    const xt::xarray<double> differentiated = xt::linalg::dot(ref.d(), ref.v());
    EXPECT_FLOATS_NEARLY_EQ(differentiated, ref.gv(), 1e-10);
  }
}
