// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oiseau/dg/field.hpp"

namespace oiseau::dg {

/// Derivatives currently represented by the first-order equation DSL.
enum class Derivative { Time, X, Y, Z };

struct DifferentialTerm {
  Field field;
  Derivative derivative = Derivative::Time;
};

struct LinearTerm {
  Field field;
  Derivative derivative = Derivative::Time;
  double coefficient = 0.0;
};

/// A linear combination of scalar fields and their first derivatives.
class LinearExpression {
 public:
  LinearExpression() = default;
  explicit LinearExpression(DifferentialTerm term)
      : m_terms{{std::move(term.field), term.derivative, 1.0}} {}

  const std::vector<LinearTerm>& terms() const { return m_terms; }

 private:
  explicit LinearExpression(std::vector<LinearTerm> terms) : m_terms(std::move(terms)) {}

  std::vector<LinearTerm> m_terms;

  friend LinearExpression operator+(const LinearExpression& lhs, const LinearExpression& rhs);
  friend LinearExpression operator+(LinearExpression&& lhs, const LinearExpression& rhs);
  friend LinearExpression operator-(LinearExpression lhs, const LinearExpression& rhs);
  friend LinearExpression operator-(LinearExpression expression);
  friend LinearExpression operator*(double coefficient, LinearExpression expression);
};

inline DifferentialTerm dt(const Field& field) { return {field, Derivative::Time}; }
inline DifferentialTerm dx(const Field& field) { return {field, Derivative::X}; }
inline DifferentialTerm dy(const Field& field) { return {field, Derivative::Y}; }
inline DifferentialTerm dz(const Field& field) { return {field, Derivative::Z}; }

inline LinearExpression operator+(const LinearExpression& lhs, const LinearExpression& rhs) {
  auto terms = lhs.terms();
  terms.insert(terms.end(), rhs.terms().begin(), rhs.terms().end());
  return LinearExpression(std::move(terms));
}

inline LinearExpression operator+(LinearExpression&& lhs, const LinearExpression& rhs) {
  if (&lhs == &rhs) {
    const auto terms = rhs.terms();
    lhs.m_terms.insert(lhs.m_terms.end(), terms.begin(), terms.end());
    return std::move(lhs);
  }
  lhs.m_terms.insert(lhs.m_terms.end(), rhs.terms().begin(), rhs.terms().end());
  return std::move(lhs);
}

inline LinearExpression operator-(LinearExpression lhs, const LinearExpression& rhs) {
  auto terms = rhs.terms();
  for (auto& term : terms) term.coefficient = -term.coefficient;
  return std::move(lhs) + LinearExpression(std::move(terms));
}

inline LinearExpression operator-(LinearExpression expression) {
  for (auto& term : expression.m_terms) term.coefficient = -term.coefficient;
  return expression;
}

inline LinearExpression operator-(DifferentialTerm term) {
  return -LinearExpression(std::move(term));
}

inline LinearExpression operator*(double coefficient, LinearExpression expression) {
  if (!std::isfinite(coefficient)) {
    throw std::invalid_argument("Equation coefficient must be finite");
  }
  for (auto& term : expression.m_terms) term.coefficient *= coefficient;
  return expression;
}

inline LinearExpression operator*(LinearExpression expression, double coefficient) {
  return coefficient * std::move(expression);
}

inline LinearExpression operator*(double coefficient, DifferentialTerm term) {
  return coefficient * LinearExpression(std::move(term));
}

inline LinearExpression operator*(DifferentialTerm term, double coefficient) {
  return coefficient * LinearExpression(std::move(term));
}

inline LinearExpression operator+(DifferentialTerm lhs, DifferentialTerm rhs) {
  return LinearExpression(std::move(lhs)) + LinearExpression(std::move(rhs));
}

inline LinearExpression operator-(DifferentialTerm lhs, DifferentialTerm rhs) {
  return LinearExpression(std::move(lhs)) - LinearExpression(std::move(rhs));
}

inline LinearExpression operator+(const LinearExpression& expression, DifferentialTerm term) {
  return expression + LinearExpression(std::move(term));
}

inline LinearExpression operator+(DifferentialTerm term, const LinearExpression& expression) {
  return LinearExpression(std::move(term)) + expression;
}

inline LinearExpression operator-(const LinearExpression& expression, DifferentialTerm term) {
  return expression - LinearExpression(std::move(term));
}

inline LinearExpression operator-(DifferentialTerm term, const LinearExpression& expression) {
  return LinearExpression(std::move(term)) - expression;
}

struct Equation {
  /// Terms normalized as `left-hand side - right-hand side == 0`.
  std::vector<LinearTerm> terms;
};

inline Equation operator==(const LinearExpression& lhs, const LinearExpression& rhs) {
  auto terms = lhs.terms();
  for (auto term : rhs.terms()) {
    term.coefficient = -term.coefficient;
    terms.push_back(std::move(term));
  }
  return {std::move(terms)};
}

inline Equation operator==(DifferentialTerm lhs, DifferentialTerm rhs) {
  return LinearExpression(std::move(lhs)) == LinearExpression(std::move(rhs));
}

inline Equation operator==(const LinearExpression& lhs, DifferentialTerm rhs) {
  return lhs == LinearExpression(std::move(rhs));
}

inline Equation operator==(DifferentialTerm lhs, const LinearExpression& rhs) {
  return LinearExpression(std::move(lhs)) == rhs;
}

/// A linear, first-order evolution system with one time derivative per field.
///
/// Each equation is normalized internally as `lhs - rhs == 0`. For a
/// conservative system `du/dt + div(F) = 0`, `flux_coefficient` returns the
/// coefficient of a field in the corresponding directional flux.
///
/// For example:
/// ```cpp
/// Field Ez{"Ez"}, Hx{"Hx"}, Hy{"Hy"};
/// EquationSystem tm({Ez, Hx, Hy}, {
///     epsilon * dt(Ez) == dx(Hy) - dy(Hx),
///     mu * dt(Hx) == -dy(Ez),
///     mu * dt(Hy) == dx(Ez),
/// });
/// ```
class EquationSystem {
 public:
  EquationSystem(std::initializer_list<Field> fields, std::initializer_list<Equation> equations)
      : EquationSystem(std::vector<Field>(fields), std::vector<Equation>(equations)) {}

  EquationSystem(std::vector<Field> fields, std::vector<Equation> equations)
      : m_fields(std::move(fields)), m_equations(std::move(equations)) {
    validate();
  }

  std::span<const Field> fields() const { return m_fields; }
  std::span<const Equation> equations() const { return m_equations; }

  /// Position of a field in states whose last axis follows `fields()` order.
  std::size_t field_index(const Field& field) const {
    const auto found = std::ranges::find_if(
        m_fields, [&field](const Field& candidate) { return candidate.same_as(field); });
    if (found != m_fields.end()) return static_cast<std::size_t>(found - m_fields.begin());
    throw std::invalid_argument("Equation references field '" + field.name() +
                                "' that is not declared in the system");
  }

  /// Coefficient in the conservative flux for `derivative` in this equation.
  double flux_coefficient(const Field& equation_field, Derivative derivative,
                          const Field& flux_field) const {
    if (derivative == Derivative::Time) {
      throw std::invalid_argument("A time derivative is not a spatial flux direction");
    }
    const auto equation_index = field_index(equation_field);
    (void)field_index(flux_field);
    const auto& equation = m_equations[m_equation_for_field[equation_index]];
    double time_coefficient = 0.0;
    double spatial_coefficient = 0.0;
    for (const auto& term : equation.terms) {
      if (term.derivative == Derivative::Time && term.field.same_as(equation_field)) {
        time_coefficient += term.coefficient;
      } else if (term.derivative == derivative && term.field.same_as(flux_field)) {
        spatial_coefficient += term.coefficient;
      }
    }
    return spatial_coefficient / time_coefficient;
  }

 private:
  void validate() {
    if (m_fields.empty() || m_fields.size() != m_equations.size()) {
      throw std::invalid_argument("EquationSystem requires one equation for each declared field");
    }
    const std::span<const Field> declared_fields(m_fields);
    for (std::size_t i = 0; i < m_fields.size(); ++i) {
      if (std::ranges::any_of(declared_fields.subspan(i + 1), [this, i](const Field& candidate) {
            return m_fields[i].same_as(candidate);
          })) {
        throw std::invalid_argument("EquationSystem contains a field more than once");
      }
    }

    std::vector<bool> has_evolution_equation(m_fields.size(), false);
    m_equation_for_field.assign(m_fields.size(), 0);
    for (std::size_t equation_index = 0; equation_index < m_equations.size(); ++equation_index) {
      const auto& equation = m_equations[equation_index];
      std::size_t time_term_count = 0;
      const Field* target = nullptr;
      for (const auto& term : equation.terms) {
        if (!std::isfinite(term.coefficient)) {
          throw std::invalid_argument("Equation coefficient must be finite");
        }
        (void)field_index(term.field);
        if (term.derivative == Derivative::Time && term.coefficient != 0.0) {
          ++time_term_count;
          target = &term.field;
        }
      }
      if (time_term_count != 1 || target == nullptr) {
        throw std::invalid_argument(
            "Each equation must contain exactly one nonzero time derivative term");
      }
      const auto target_index = field_index(*target);
      if (has_evolution_equation[target_index]) {
        throw std::invalid_argument("EquationSystem has multiple evolution equations for field '" +
                                    target->name() + "'");
      }
      has_evolution_equation[target_index] = true;
      m_equation_for_field[target_index] = equation_index;
    }
    if (std::ranges::find(has_evolution_equation, false) != has_evolution_equation.end()) {
      throw std::invalid_argument("EquationSystem is missing an evolution equation for a field");
    }
  }

  std::vector<Field> m_fields;
  std::vector<Equation> m_equations;
  std::vector<std::size_t> m_equation_for_field;
};

}  // namespace oiseau::dg
