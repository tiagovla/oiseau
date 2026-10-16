// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"

/**
 * @file ref_tetrahedron.hpp
 * @brief Defines the reference tetrahedron element used in nodal Discontinuous Galerkin methods.
 *
 * This file contains the declaration of the `RefTetrahedron` class, which represents
 * a reference tetrahedron element, and associated internal helper functions in the `detail`
 * namespace for generating nodes and performing coordinate transformations.
 */

namespace oiseau::dg::nodal {

/**
 * @class RefTetrahedron
 * @brief Represents a reference tetrahedron element for nodal Discontinuous Galerkin methods.
 *
 * The reference tetrahedron is defined in a canonical (r, s, t) coordinate system,
 * typically mapping to the region { (r,s,t) | r,s,t >= -1, r+s+t <= -1 }.
 * This class computes properties like nodal distribution, Vandermonde matrices,
 * and differentiation operators for a specified polynomial order.
 */
class RefTetrahedron : public RefElement {
 public:
  /**
   * @brief Constructs a RefTetrahedron object of a specific polynomial order.
   * @param order The polynomial order for the basis functions on the tetrahedron.
   *              Determines the number of nodes and approximation accuracy.
   *              The total degree of basis polynomials will be at most `order`.
   */
  explicit RefTetrahedron(unsigned order);

  /**
   * @brief Evaluates a 3D basis function on the reference tetrahedron.
   *
   * The basis functions are typically orthogonal polynomials (e.g., a product of
   * Jacobi polynomials) defined in a transformed coordinate system (a,b,c)
   * derived from the reference (r,s,t) coordinates. The specific form and
   * the range of indices i,j,k may depend on `this->order()`.
   * Typically, \f$i+j+k \le \text{order}\f$.
   *
   * @param rst 2D array (shape: N_points × 3) of (r, s, t) coordinates where to evaluate.
   * @param i Polynomial degree index (e.g., for the first transformed coordinate).
   * @param j Polynomial degree index (e.g., for the second transformed coordinate).
   * @param k Polynomial degree index (e.g., for the third transformed coordinate).
   * @return  1D array (length N_points) of basis function values at each (r, s, t) point.
   */
  static xt::xarray<double> basis_function(const xt::xarray<double> &rst, int i, int j, int k);

  /**
   * @brief Computes the gradient of a 3D basis function on the reference tetrahedron.
   *
   * Computes derivatives of the basis function \f$\phi_{ijk}\f$ with respect to
   * the reference coordinates (r, s, t).
   *
   * @param rst 2D array (shape: N_points × 3) of (r, s, t) coordinates.
   * @param i Polynomial degree index.
   * @param j Polynomial degree index.
   * @param k Polynomial degree index.
   * @return  2D array (shape: N_points × 3) containing (∂/∂r, ∂/∂s, ∂/∂t) for each point.
   */
  static xt::xarray<double> grad_basis_function(const xt::xarray<double> &rst, int i, int j, int k);

 private:
  /**
   * @brief Computes the Vandermonde matrix for the reference tetrahedron.
   *
   * The Vandermonde matrix V has entries V_{p,α} = φ_α(r_p, s_p, t_p), where
   * (r_p, s_p, t_p) are the nodal points and φ_α are the basis functions
   * indexed by a multi-index α=(i,j,k).
   *
   * @param rst 2D array (shape: N_nodes × 3) of (r, s, t) nodal coordinates.
   * @return    2D array (shape: N_nodes × N_basis) representing the Vandermonde matrix.
   */
  xt::xarray<double> vandermonde(const xt::xarray<double> &rst) const override;

  /**
   * @brief Computes the gradient Vandermonde matrix for the reference tetrahedron.
   *
   * The gradient Vandermonde matrix G_V has entries:
   * G_{p,α,0} = ∂_r φ_α(r_p,s_p,t_p)
   * G_{p,α,1} = ∂_s φ_α(r_p,s_p,t_p)
   * G_{p,α,2} = ∂_t φ_α(r_p,s_p,t_p)
   *
   * @param rst 2D array (shape: N_nodes × 3) of (r, s, t) nodal coordinates.
   * @return    3D array (N_nodes × N_basis × 3) of basis function gradients at nodal points.
   */
  xt::xarray<double> grad_vandermonde(const xt::xarray<double> &rst) const;

  /**
   * @brief Computes the differentiation operator matrices Dr, Ds, Dt.
   *
   * Dr = V_r * V^{-1}, Ds = V_s * V^{-1}, Dt = V_t * V^{-1} at the nodal points,
   * where V_r, V_s, V_t are slices of the gradient Vandermonde matrix.
   *
   * @param v  2D array (N_nodes × N_basis) Vandermonde matrix.
   * @param gv 3D array (N_nodes × N_basis × 3) gradient Vandermonde matrix.
   * @return   3D array (N_nodes × N_basis × 3) containing Dr ([:,:,0]), Ds ([:,:,1]), and Dt
   * ([:,:,2]). Note: The output shape is typically N_nodes x N_nodes x 3 if V is square (N_nodes =
   * N_basis). Or, if it's N_nodes x N_basis x 3, it implies operators that map modal coefficients
   * to nodal derivatives. Assuming N_nodes = N_basis, returns differentiation matrices of shape
   * (N_nodes x N_nodes x 3).
   */
  static xt::xarray<double> grad_operator(const xt::xarray<double> &v,
                                          const xt::xarray<double> &gv);
};

namespace detail {

/**
 * @brief Transforms coordinates from a canonical "equilateral" tetrahedron in XYZ space
 *        to the reference (r,s,t) computational space.
 *
 * The "equilateral" tetrahedron might be one with vertices like (0,0,0), (1,0,0), (0.5, sqrt(3)/2,
 * 0), (0.5, sqrt(3)/6, sqrt(6)/3). The reference (r,s,t) space is typically { (r,s,t) | r,s,t >=
 * -1, r+s+t <= -1 }.
 *
 * @param coords 2D array (shape: N_points × 3) of XYZ coordinates.
 * @return       2D array (shape: N_points × 3) of (r,s,t) coordinates.
 */
xt::xarray<double> equilateral_xyz_to_rst(const xt::xarray<double> &coords);

/**
 * @brief Transforms coordinates from the reference (r,s,t) space to a auxiliary (a,b,c) space.
 *
 * The (a,b,c) coordinates are often collapsed coordinates (e.g., Duffy-like or Koornwinder-type)
 * in which the basis functions (e.g., products of Jacobi polynomials) are conveniently defined.
 * For example, for a tetrahedron, these might be:
 * \f$a = r\f$
 * \f$b = \frac{2(1+s)}{1-t} - 1\f$ (if t != 1)
 * \f$c = t\f$
 * (This specific example is for a pyramid, tetrahedral transformations are similar).
 *
 * @param rst 2D array (shape: N_points × 3) of (r,s,t) coordinates.
 * @return    2D array (shape: N_points × 3) of (a,b,c) coordinates.
 */
xt::xarray<double> rst_to_abc(const xt::xarray<double> &rst);

/**
 * @brief Generates nodal points on the reference tetrahedron for a given polynomial order.
 *
 * These nodes are typically optimized for polynomial interpolation, often generated
 * using methods like Warp & Blend (Warburton, 2006).
 *
 * @param p The polynomial order. The number of nodes will be (p+1)(p+2)(p+3)/6.
 * @return  2D array ( (p+1)(p+2)(p+3)/6 × 3 ) of (r,s,t) node coordinates.
 */
xt::xarray<double> generate_tetrahedron_nodes(unsigned p);

/**
 * @brief Generates approximately equidistant nodes on the reference tetrahedron.
 *
 * These nodes are typically used as a starting point for warping algorithms (like Warp & Blend)
 * or for low-order approximations.
 *
 * @param n The polynomial order, determining the density of points.
 *          The number of points is (n+1)(n+2)(n+3)/6.
 * @return  2D array ( (n+1)(n+2)(n+3)/6 × 3 ) of (r,s,t) node coordinates.
 */
xt::xarray<double> generate_tetrahedron_equidistant_nodes(unsigned n);

/**
 * @brief Helper function for the Warp & Blend algorithm: computes node shifts
 *        for points near the faces of the tetrahedron.
 *
 * This function is part of the process to improve the quality of nodal distributions
 * on the tetrahedron by "warping" an initial set of points.
 *
 * @param p Polynomial order.
 * @param pval A parameter, possibly related to the evaluation point's properties or a scaling
 * factor.
 * @param l2 Barycentric coordinate \f$\zeta_2\f$.
 * @param l3 Barycentric coordinate \f$\zeta_3\f$.
 * @param l4 Barycentric coordinate \f$\zeta_4\f$.
 * @return   Shift values (xt::xarray<double>) for the nodes.
 */
xt::xarray<double> warp_shift_face_3d(int p, double pval, const xt::xarray<double> &l2,
                                      const xt::xarray<double> &l3, const xt::xarray<double> &l4);
/**
 * @brief Helper function for the Warp & Blend algorithm: evaluates a shifting polynomial.
 *
 * This polynomial contributes to the overall warp transformation, often based on distances
 * to element boundaries (represented by barycentric coordinates).
 *
 * @param p Polynomial order.
 * @param pval A parameter, possibly related to the evaluation point or a scaling factor for the
 * shift.
 * @param l1 Barycentric coordinate \f$\zeta_1\f$.
 * @param l2 Barycentric coordinate \f$\zeta_2\f$.
 * @param l3 Barycentric coordinate \f$\zeta_3\f$.
 * @return   Evaluated shift values (xt::xarray<double>).
 */
xt::xarray<double> eval_shift(int p, double pval, const xt::xarray<double> &l1,
                              const xt::xarray<double> &l2, const xt::xarray<double> &l3);
/**
 * @brief Helper function for the Warp & Blend algorithm: evaluates the warp factor or displacement.
 *
 * This function computes the actual warp displacement for a set of points `xout`
 * based on the polynomial order `p` and potentially a reference set of 1D nodes `xnodes`
 * (e.g., alpha-optimized Gauss-Lobatto nodes).
 *
 * @param p Polynomial order.
 * @param xnodes 1D array of reference nodes (e.g., roots of a Jacobi polynomial) used to define the
 * warp.
 * @param xout   1D array of points (e.g., barycentric coordinates) at which to evaluate the warp
 * factor.
 * @return       Warp factors or displacements (xt::xarray<double>).
 */
xt::xarray<double> eval_warp(int p, const xt::xarray<double> &xnodes,
                             const xt::xarray<double> &xout);

}  // namespace detail

}  // namespace oiseau::dg::nodal
