// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <stdexcept>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/test_macros.hpp"
#include "oiseau/utils/integration.hpp"

namespace oiseau {
namespace {

namespace integration = utils::integration;

TEST(IntegrationTest, GaussQuadratureIntegratesMonomials) {
  constexpr int order = 50;
  const auto [nodes, weights] = integration::quadrature(order);
  EXPECT_EQ(nodes.size(), 26);
  for (int degree = 0; degree <= order; ++degree) {
    const double expected = degree % 2 == 0 ? 2.0 / (degree + 1) : 0.0;
    const double integral = xt::sum(xt::pow(nodes, degree) * weights)();
    EXPECT_NEAR(integral, expected, 1e-14) << "Degree " << degree;
  }
}

TEST(IntegrationTest, TriangleCubatureIntegratesMonomials) {
  constexpr int order = 40;
  const auto [nodes, weights] = integration::triangle_cubature(order);
  const auto x = xt::col(nodes, 0);
  const auto y = xt::col(nodes, 1);

  EXPECT_NEAR(xt::sum(weights)(), 2.0, 5e-14);
  for (int degree = 1; degree <= order; ++degree) {
    const double expected = degree % 2 == 0 ? 2.0 / (degree + 1) : -2.0 / (degree + 2);
    const double x_integral = xt::sum(xt::pow(x, degree) * weights)();
    const double y_integral = xt::sum(xt::pow(y, degree) * weights)();
    EXPECT_NEAR(x_integral, expected, 5e-14) << "x degree " << degree;
    EXPECT_NEAR(y_integral, expected, 5e-14) << "y degree " << degree;
  }
}

TEST(IntegrationTest, TriangleCubatureRejectsNonPositiveOrders) {
  EXPECT_THROW(integration::triangle_cubature(-1), std::out_of_range);
  EXPECT_THROW(integration::triangle_cubature(0), std::out_of_range);
}

TEST(IntegrationTest, QuadratureRejectsNonPositiveOrders) {
  EXPECT_THROW(integration::quadrature(-1), std::out_of_range);
  EXPECT_THROW(integration::quadrature(0), std::out_of_range);
}

TEST(IntegrationTest, GaussJacobiNodesAndWeights) {
  const auto [nodes, weights] = integration::jacobi_gq(3, 1.0, 2.0);
  const xt::xarray<double> expected_nodes = {-0.65077886, -0.15637043, 0.37348938, 0.79729627};
  const xt::xarray<double> expected_weights = {0.08666291, 0.44123335, 0.59015336, 0.21528371};
  EXPECT_FLOATS_NEARLY_EQ(expected_nodes, nodes, 1e-4);
  EXPECT_FLOATS_NEARLY_EQ(expected_weights, weights, 1e-4);
}

TEST(IntegrationTest, OrderOneLobattoJacobiNodesAreEndpoints) {
  const auto nodes = integration::jacobi_gl(1, 0.0, 0.0);
  const xt::xarray<double> expected_nodes = {-1.0, 1.0};
  EXPECT_FLOATS_NEARLY_EQ(expected_nodes, nodes, 1e-10);
}

TEST(IntegrationTest, LobattoJacobiInteriorNodes) {
  const auto nodes = integration::jacobi_gl(4, 1.0, 2.0);
  const xt::xarray<double> expected_nodes = {-1.0, -0.4611587, 0.10438533, 0.62950065, 1.0};
  EXPECT_FLOATS_NEARLY_EQ(expected_nodes, nodes, 1e-4);
}

}  // namespace
}  // namespace oiseau
