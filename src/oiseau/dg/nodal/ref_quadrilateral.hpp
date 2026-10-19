// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"

/**
 * @file ref_quadrilateral.hpp
 * @brief Defines the reference quadrilateral element used in nodal Discontinuous Galerkin methods.
 *
 * This file contains the declaration of the `RefQuadrilateral` class, which represents
 * a reference square element, and associated internal helper functions in the `detail`
 * namespace for generating nodes and performing coordinate transformations.
 *
 * Node layout (order = 2):
 *
 *   s ↑
 *     |
 *     | (6)——(7)——(8)
 *     |  |     |     |
 *     | (3)——(4)——(5)
 *     |  |     |     |
 *     | (0)——(1)——(2)
 *     +——————→ r
 *
 *   Nodes are numbered in lexicographic order from bottom-left to top-right.
 */

namespace oiseau::dg::nodal {

/**
 * @class RefQuadrilateral
 * @brief Represents a reference quadrilateral element for nodal Discontinuous Galerkin methods.
 *
 * The reference quadrilateral is defined in (r, s) coordinates on the square
 * [-1, 1] × [-1, 1]. This class computes properties like nodal distribution,
 * Vandermonde matrices, and differentiation operators for a specified polynomial order.
 */
class RefQuadrilateral : public RefElement {
 public:
  /**
   * @brief Constructs a RefQuadrilateral object of a specific polynomial order.
   * @param order The polynomial order for the basis functions on the quadrilateral.
   *              Determines the number of nodes in each direction and approximation accuracy.
   */
  explicit RefQuadrilateral(unsigned order);

  /**
   * @brief Static member function that evaluates a 2D tensor-product basis function
   *        on the reference quadrilateral.
   *
   * The basis function is defined as the product of 1D Jacobi polynomials:
   * \f$\phi_{ij}(r, s) = P_i^{0,0}(r) \; P_j^{0,0}(s)\f$,
   * evaluated at (r,s) points.
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates where to evaluate.
   * @param i  Polynomial degree index in the r-direction (0 ≤ i ≤ order).
   * @param j  Polynomial degree index in the s-direction (0 ≤ j ≤ order).
   * @return   1D array (length N_points) of basis function values at each (r, s) point.
   */
  static xt::xarray<double> basis_function(const xt::xarray<double>& rs, int i, int j);

  /**
   * @brief Static member function that computes the gradient of a tensor-product
   *        basis function on the reference quadrilateral.
   *
   * Computes derivatives with respect to r and s separately using:
   * \f$\partial_r \phi_{ij} = P'_i(r) P_j(s),\quad
   *   \partial_s \phi_{ij} = P_i(r) P'_j(s)\f$.
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
   * @param i  Polynomial degree index in the r-direction.
   * @param j  Polynomial degree index in the s-direction.
   * @return   2D array (shape: N_points × 2) containing (∂/∂r, ∂/∂s) for each point.
   */
  static xt::xarray<double> grad_basis_function(const xt::xarray<double>& rs, int i, int j);

 private:
  /**
   * @brief Computes the Vandermonde matrix for the reference quadrilateral.
   *
   * V has entries V_{k,α} = φ_α(r_k, s_k), where (r_k, s_k) are the nodal points
   * and φ_α are the tensor-product basis functions indexed by multi-index α.
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
   * @return   2D array (shape: N_points × N_basis) representing the Vandermonde matrix.
   */
  xt::xarray<double> vandermonde(const xt::xarray<double>& rs) const override;

  /**
   * @brief Computes the gradient Vandermonde matrix for the reference quadrilateral.
   *
   * G_V has entries G_{k,α,0} = ∂_r φ_α(r_k,s_k), G_{k,α,1} = ∂_s φ_α(r_k,s_k).
   *
   * @param rs 2D array (shape: N_points × 2) of (r, s) coordinates.
   * @return   3D array (N_points × N_basis × 2) of gradients.
   */
  xt::xarray<double> grad_vandermonde(const xt::xarray<double>& rs) const;

  /**
   * @brief Computes the differentiation operator matrices Dr and Ds.
   *
   * Dr = V_r * V^{-1}, Ds = V_s * V^{-1} at the nodal points,
   * where V_r and V_s are the derivative Vandermonde slices.
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
 * @brief Generates tensor-product Gauss-Lobatto nodes on the reference quadrilateral.
 *
 * Computes the 1D Gauss-Lobatto points of the given order and forms the tensor grid.
 *
 * @param order Polynomial order (number of points per direction = order+1).
 * @return      2D array ( (order+1)^2 × 2 ) of (r, s) node coordinates.
 */
xt::xarray<double> generate_quadrilateral_nodes(unsigned order);

}  // namespace detail
}  // namespace oiseau::dg::nodal
