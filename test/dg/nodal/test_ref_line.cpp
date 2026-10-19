// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_line.hpp"
#include "oiseau/dg/nodal/utils.hpp"
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

TEST(test_ref_line, vandermonde_interpolation) {
  for (unsigned order = 3; order <= 5; ++order) {
    // Coarse and fine reference elements
    auto coarse_ref = oiseau::dg::nodal::RefLine(order);
    auto fine_ref = oiseau::dg::nodal::RefLine(order + 2);

    // Vandermonde matrices at coarse and fine nodes
    auto V_coarse = coarse_ref.v();
    auto V_fine = oiseau::dg::nodal::utils::vandermonde_1d(order, fine_ref.r());

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
