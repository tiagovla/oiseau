// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"

namespace oiseau::dg::nodal {

/**
 * @class RefLine
 * @brief Represents a reference line element used in nodal Discontinuous Galerkin methods.
 *
 * This class implements the necessary basis functions, their gradients, and supporting
 * operators for a 1D reference element (line) of arbitrary polynomial order.
 *
 * ASCII representation of the reference line node numbering (order = 4):
 *
 * (0)   (1)   (2)   (3)   (4)
 *  |-----|-----|-----|-----|
 * -1                        1
 *
 * The line spans the reference interval [-1, 1]. Nodes are placed at Gauss, Gauss–Lobatto,
 * or Jacobi–Gauss points depending on the implementation in `generate_line_nodes`.
 */
class RefLine : public RefElement {
 public:
  /**
   * @brief Constructs a RefLine element of the given polynomial order.
   * @param order The polynomial order of the reference element.
   */
  explicit RefLine(unsigned order);

  /**
   * @brief Evaluates the i-th Lagrange basis function at the given reference points.
   * @param r Reference coordinates where the basis function is evaluated.
   * @param i Index of the basis function to evaluate.
   * @return The value of the i-th basis function at each point in r.
   */
  static xt::xarray<double> basis_function(const xt::xarray<double>& r, int i);

  /**
   * @brief Evaluates the gradient of the i-th Lagrange basis function at the given reference
   * points.
   * @param r Reference coordinates where the gradient is evaluated.
   * @param i Index of the basis function.
   * @return The gradient of the i-th basis function at each point in r.
   */
  static xt::xarray<double> grad_basis_function(const xt::xarray<double>& r, int i);

 private:
  /**
   * @brief Computes the Vandermonde matrix for the given reference coordinates.
   * @param r Reference coordinates.
   * @return The Vandermonde matrix used for interpolation.
   */
  xt::xarray<double> vandermonde(const xt::xarray<double>& r) const override;

  /**
   * @brief Computes the gradient of the Vandermonde matrix for the given reference coordinates.
   * @param r Reference coordinates.
   * @return The gradient Vandermonde matrix.
   */
  xt::xarray<double> grad_vandermonde(const xt::xarray<double>& r) const;

  /**
   * @brief Computes the differentiation matrix from Vandermonde and its gradient.
   * @param v The Vandermonde matrix.
   * @param gv The gradient of the Vandermonde matrix.
   * @return The gradient (differentiation) operator matrix.
   */
  static xt::xarray<double> grad_operator(const xt::xarray<double>& v,
                                          const xt::xarray<double>& gv);
};

namespace detail {

/**
 * @brief Generates nodal points on the reference line element for a given polynomial order.
 *
 * The generated nodes typically lie in the interval [-1, 1] and may correspond to
 * Gauss, Gauss–Lobatto, or Jacobi–Gauss quadrature points depending on the implementation.
 *
 * @param order Polynomial order of the element.
 * @return An array of nodal coordinates on the reference line.
 */
xt::xarray<double> generate_line_nodes(unsigned order);

}  // namespace detail
}  // namespace oiseau::dg::nodal
