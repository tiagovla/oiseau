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

#include "oiseau/dg/nodal/ref_tetrahedron.hpp"
#include "oiseau/dg/nodal/utils.hpp"
#include "oiseau/test_macros.hpp"

TEST(test_ref_tetrahedron, get_order) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefTetrahedron(order);
    EXPECT_EQ(ref.order(), order);
  }
}

TEST(test_ref_tetrahedron, get_number_of_nodes) {
  for (unsigned order = 1; order <= 5; ++order) {
    auto ref = oiseau::dg::nodal::RefTetrahedron(order);
    unsigned n = (order + 1) * (order + 2) * (order + 3) / 6;
    EXPECT_EQ(ref.number_of_nodes(), n);
  }
}

TEST(test_ref_tetrahedron, vandermonde_interpolation) {
  for (unsigned order = 3; order <= 5; ++order) {
    // Coarse and fine reference elements
    auto coarse_ref = oiseau::dg::nodal::RefTetrahedron(order);
    auto fine_ref = oiseau::dg::nodal::RefTetrahedron(order + 2);

    // Vandermonde matrices at coarse and fine nodes
    auto V_coarse = coarse_ref.v();
    auto V_fine = oiseau::dg::nodal::utils::vandermonde_3d(order, fine_ref.r());

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
