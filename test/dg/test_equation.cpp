// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "oiseau/dg/equation.hpp"

namespace {

using oiseau::dg::Derivative;
using oiseau::dg::dt;
using oiseau::dg::dx;
using oiseau::dg::dy;
using oiseau::dg::dz;
using oiseau::dg::Equation;
using oiseau::dg::EquationSystem;
using oiseau::dg::Field;

}  // namespace

TEST(test_dg_equation, tm_maxwell_flux_matrices_follow_declared_field_order) {
  constexpr double epsilon = 2.0;
  constexpr double mu = 4.0;
  Field Ez{"Ez"}, Hx{"Hx"}, Hy{"Hy"};
  EquationSystem equations({Ez, Hx, Hy}, {
                                             epsilon * dt(Ez) == dx(Hy) - dy(Hx),
                                             mu * dt(Hx) == -dy(Ez),
                                             mu * dt(Hy) == dx(Ez),
                                         });

  const std::array<Field, 3> fields = {Ez, Hx, Hy};
  const std::array<std::array<double, 3>, 3> expected_x = {
      {{0.0, 0.0, -0.5}, {0.0, 0.0, 0.0}, {-0.25, 0.0, 0.0}}};
  const std::array<std::array<double, 3>, 3> expected_y = {
      {{0.0, 0.5, 0.0}, {0.25, 0.0, 0.0}, {0.0, 0.0, 0.0}}};

  EXPECT_EQ(equations.field_index(Ez), 0);
  EXPECT_EQ(equations.field_index(Hx), 1);
  EXPECT_EQ(equations.field_index(Hy), 2);
  for (std::size_t row = 0; row < fields.size(); ++row) {
    for (std::size_t column = 0; column < fields.size(); ++column) {
      EXPECT_DOUBLE_EQ(equations.flux_coefficient(fields[row], Derivative::X, fields[column]),
                       expected_x[row][column]);
      EXPECT_DOUBLE_EQ(equations.flux_coefficient(fields[row], Derivative::Y, fields[column]),
                       expected_y[row][column]);
    }
  }
}

TEST(test_dg_equation, maps_equations_to_fields_when_equation_order_differs) {
  Field u{"u"}, v{"v"};
  EquationSystem equations({u, v}, {dt(v) == dz(u), dt(u) == dx(v)});

  EXPECT_EQ(equations.field_index(u), 0);
  EXPECT_EQ(equations.field_index(v), 1);
  EXPECT_DOUBLE_EQ(equations.flux_coefficient(u, Derivative::X, v), -1.0);
  EXPECT_DOUBLE_EQ(equations.flux_coefficient(v, Derivative::Z, u), -1.0);
}

TEST(test_dg_equation, combines_spatial_terms_and_divides_by_time_coefficient) {
  Field u{"u"}, v{"v"};
  EquationSystem equations({u, v}, {
                                       2.0 * dt(u) == 3.0 * dx(v) - dy(u) + dx(v),
                                       4.0 * dt(v) == -dy(u),
                                   });

  EXPECT_DOUBLE_EQ(equations.flux_coefficient(u, Derivative::X, v), -2.0);
  EXPECT_DOUBLE_EQ(equations.flux_coefficient(u, Derivative::Y, u), 0.5);
  EXPECT_DOUBLE_EQ(equations.flux_coefficient(v, Derivative::Y, u), 0.25);
}

TEST(test_dg_equation, rejects_empty_or_mismatched_field_and_equation_lists) {
  Field u{"u"};

  EXPECT_THROW((EquationSystem(std::vector<Field>{}, std::vector<Equation>{})),
               std::invalid_argument);
  EXPECT_THROW((EquationSystem(std::vector<Field>{u}, std::vector<Equation>{})),
               std::invalid_argument);
}

TEST(test_dg_equation, rejects_duplicate_fields_and_undeclared_references) {
  Field u{"u"}, v{"v"};

  EXPECT_THROW((EquationSystem({u, u}, {dt(u) == dx(u), dt(u) == dy(u)})), std::invalid_argument);
  EXPECT_THROW((EquationSystem({u}, {dt(u) == dx(v)})), std::invalid_argument);
}

TEST(test_dg_equation, requires_exactly_one_nonzero_time_derivative_per_equation) {
  Field u{"u"}, v{"v"};

  EXPECT_THROW((EquationSystem({u}, {dx(u) == dy(u)})), std::invalid_argument);
  EXPECT_THROW((EquationSystem({u}, {0.0 * dt(u) == dx(u)})), std::invalid_argument);
  EXPECT_THROW((EquationSystem({u, v}, {dt(u) + dt(v) == dx(u), dt(v) == dx(v)})),
               std::invalid_argument);
}

TEST(test_dg_equation, rejects_multiple_evolution_equations_for_one_field) {
  Field u{"u"}, v{"v"};

  EXPECT_THROW((EquationSystem({u, v}, {dt(u) == dx(v), dt(u) == dy(v)})), std::invalid_argument);
}

TEST(test_dg_equation, rejects_nonfinite_coefficients) {
  Field u{"u"};
  const auto infinity = std::numeric_limits<double>::infinity();

  EXPECT_THROW(infinity * dt(u), std::invalid_argument);

  Equation equation;
  equation.terms.push_back({u, Derivative::Time, infinity});
  EXPECT_THROW((EquationSystem(std::vector<Field>{u}, std::vector<Equation>{equation})),
               std::invalid_argument);
}

TEST(test_dg_equation, field_identity_and_flux_direction_are_validated) {
  Field u{"u"}, another_u_with_same_name{"u"};
  EquationSystem equations({u}, {dt(u) == dx(u)});

  EXPECT_THROW(equations.field_index(another_u_with_same_name), std::invalid_argument);
  EXPECT_THROW(equations.flux_coefficient(another_u_with_same_name, Derivative::X, u),
               std::invalid_argument);
  EXPECT_THROW(equations.flux_coefficient(u, Derivative::X, another_u_with_same_name),
               std::invalid_argument);
  EXPECT_THROW(equations.flux_coefficient(u, Derivative::Time, u), std::invalid_argument);
}
