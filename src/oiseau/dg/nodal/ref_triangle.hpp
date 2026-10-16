// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "oiseau/dg/nodal/ref_element.hpp"
#include "xtensor/core/xtensor_forward.hpp"

/**
 * @file ref_triangle.hpp
 * @brief Defines the reference triangle element used in nodal Discontinuous Galerkin methods.
 *
 * This file contains the declaration of the `RefTriangle` class, which represents
 * a reference triangular element, and associated internal helper functions in the `detail`
 * namespace for generating nodes, basis functions, and performing coordinate transformations.

 * ASCII representation of the reference triangle node numbering (order=4):
 *
 * (14)
 *  |   \
 * (12)  (13)
 *  |        \
 * (9)  (10)  (11)
 *  |            \
 * (5)  (6)  (7)  (8)
 *  |                 \
 * (0)--(1)--(2)--(3)--(4)
 */

namespace oiseau::dg::nodal {

/**
 * @class RefTriangle
 * @brief Represents a reference triangular element for nodal Discontinuous Galerkin methods.
 *
 * The reference triangle is defined in (r, s) coordinates with vertices at
 * (-1, -1), (1, -1), and (-1, 1). This class computes and stores properties
 * such as nodal distribution, Vandermonde matrix, and differentiation matrices
 * for a specified polynomial order.
 */
class RefTriangle : public RefElement {
 public:
  /**
   * @brief Constructs a RefTriangle object of a specific polynomial order.
   * @param order The polynomial order for the basis functions on the triangle.
   *              Determines the number of nodes and approximation accuracy.
   */
  explicit RefTriangle(unsigned order);

  /**
   * @brief Static member function that evaluates a 2D orthonormal basis function
   *        on the reference triangle in collapsed (a, b) coordinates.
   *
   * The basis function is defined by:
   * \f$\sqrt{2}\,P_i^{0,0}(a)\,P_j^{2i+1,0}(b)\,(1 - b)^i\f$
   * where \f$P_n^{\alpha,\beta}\f$ are Jacobi polynomials.
   *
   * @param ab 2D array (shape: N_points × 2) of (a, b) coordinates where to evaluate the basis.
   * @param i  Polynomial degree index for the a-direction component.
   * @param j  Polynomial degree index for the b-direction component.
   * @return   1D array (length N_points) of basis function values at each (a, b) point.
   */
  static xt::xarray<double> basis_function(const xt::xarray<double>& ab, int i, int j);

  /**
   * @brief Static member function that computes the gradient of a 2D orthonormal
   *        basis function on the reference triangle in (r, s) coordinates.
   *
   * Uses the chain rule to transform derivatives from collapsed (a, b) to (r, s).
   *
   * @param ab 2D array (shape: N_points × 2) of (a, b) coordinates from `rs_to_ab`.
   * @param i  Polynomial degree index for the a-component.
   * @param j  Polynomial degree index for the b-component.
   * @return   2D array (shape: N_points × 2) containing (∂/∂r, ∂/∂s) for each point.
   */
  static xt::xarray<double> grad_basis_function(const xt::xarray<double>& ab, int i, int j);

 private:
  /**
   * @brief Computes the Vandermonde matrix for the reference triangle.
   *
   * The Vandermonde matrix V has entries V_{k,j} = φ_j(r_k, s_k), where
   * (r_k, s_k) are the nodal points and φ_j are the basis functions.
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
   * @return   2D array (shape: N_points × N_basis) representing the Vandermonde matrix.
   */
  xt::xarray<double> vandermonde(const xt::xarray<double>& rs) const override;

  /**
   * @brief Computes the gradient of the Vandermonde matrix for the reference triangle.
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
   * @return   3D array (shape: N_points × N_basis × 2) of gradients (∂/∂r, ∂/∂s).
   */
  xt::xarray<double> grad_vandermonde(const xt::xarray<double>& rs) const;

  /**
   * @brief Computes the differentiation operator matrices Dr and Ds.
   *
   * Dr = V_r * V^{-1}, Ds = V_s * V^{-1}, where V_r and V_s are the parts
   * of the gradient Vandermonde matrix and V^{-1} is the inverse Vandermonde
   * at element nodes.
   *
   * @param v  2D array (N_nodes × N_basis) Vandermonde matrix.
   * @param gv 3D array (N_nodes × N_basis × 2) gradient Vandermonde.
   * @return   3D array (N_nodes × N_basis × 2) containing Dr ([:,:,0]) and Ds ([:,:,1]).
   */
  static xt::xarray<double> grad_operator(const xt::xarray<double>& v,
                                          const xt::xarray<double>& gv);
};

namespace detail {

/**
 * @brief Converts (r, s) coordinates to collapsed (a, b) coordinates for the triangle.
 *
 * a = 2*(1+r)/(1-s) - 1,  b = s; with special handling when s == 1.
 *
 * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
 * @return   2D array (shape: N_points × 2) of (a, b) coordinates.
 */
xt::xarray<double> rs_to_ab(const xt::xarray<double>& rs);

/**
 * @brief Generates nodal points on an equilateral triangle via Warp & Blend.
 *
 * The result is then mapped to the (r, s) reference triangle.
 *
 * @param order Polynomial order (number of points = (order+1)*(order+2)/2).
 * @return      2D array (N_points × 2) of (x, y) on the equilateral triangle.
 */
xt::xarray<double> generate_triangle_nodes(unsigned order);

/**
 * @brief Maps equilateral triangle (x, y) to reference (r, s) coordinates.
 *
 * @param xy 2D array (N_points × 2) of equilateral (x, y).
 * @return   2D array (N_points × 2) of reference (r, s).
 */
xt::xarray<double> equilateral_xy_to_rs(const xt::xarray<double>& xy);

/**
 * @brief Computes 1D warp factor for node distribution on triangle edges.
 *
 * @param order Polynomial order.
 * @param rs    1D array (length = order+1) of coordinates along an edge.
 * @return      1D array (length = order+1) of warp offsets.
 */
xt::xarray<double> warp_factor(unsigned order, const xt::xarray<double>& rs);

}  // namespace detail
}  // namespace oiseau::dg::nodal
