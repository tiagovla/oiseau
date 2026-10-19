// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>

namespace oiseau::dg::nodal::utils {

/**
 * @brief calculates the gauss-lobatto quadrature points and weights
 *
 * @param n order of the quadrature
 * @param alpha parameter of the jacobi polynomial
 * @param beta parameter of the jacobi polynomial
 */
std::pair<xt::xarray<double>, xt::xarray<double>> jacobi_gq(unsigned n, double alpha, double beta);

/**
 * @brief calculates the gauss-legendre quadrature points
 *
 * @param n order of the quadrature
 * @param alpha parameter of the jacobi polynomial
 * @param beta parameter of the jacobi polynomial
 */
xt::xarray<double> jacobi_gl(unsigned n, double alpha, double beta);

xt::xarray<double> d_matrix_1d(const xt::xarray<double> &v, const xt::xarray<double> &gv);
xt::xarray<double> d_matrix_1d(unsigned n, const xt::xarray<double> &r);

/**
 * @brief initialize the gradient of the modal basis at r at order n
 *
 * @param n order
 * @param r locations [-1, 1]
 */
xt::xarray<double> vandermonde_1d(unsigned n, const xt::xarray<double> &r);
xt::xarray<double> grad_vandermonde_1d(unsigned n, const xt::xarray<double> &r);

xt::xarray<double> d_matrix_2d(const xt::xarray<double> &v, const xt::xarray<double> &gv);
xt::xarray<double> vandermonde_2d(unsigned n, const xt::xarray<double> &rs);
xt::xarray<double> vandermonde_2d_tensor(unsigned n, const xt::xarray<double> &rs);
xt::xarray<double> grad_vandermonde_2d(unsigned n, const xt::xarray<double> &rs);
xt::xarray<double> grad_vandermonde_2d_tensor(unsigned n, const xt::xarray<double> &rs);

xt::xarray<double> simplexp_2d(const xt::xarray<double> &ab, int i, int j);
xt::xarray<double> tensorp_2d(const xt::xarray<double> &rs, int i, int j);
xt::xarray<double> grad_simplexp_2d(const xt::xarray<double> &ab, int i, int j);
xt::xarray<double> grad_tensorp_2d(const xt::xarray<double> &rs, int i, int j);

xt::xarray<double> generate_triangle_equidistant_nodes(unsigned n);
xt::xarray<double> warp_factor(unsigned n, const xt::xarray<double> &rout);
xt::xarray<double> conversion_equilateral_xy_to_rs(const xt::xarray<double> &coords);
xt::xarray<double> conversion_rs_to_ab(const xt::xarray<double> &rs);
xt::xarray<double> generate_triangle_nodes(unsigned n);
xt::xarray<double> generate_tetrahedron_equidistant_nodes(unsigned n);
xt::xarray<double> generate_quadrilateral_nodes(unsigned n);
xt::xarray<double> generate_hexahedron_nodes(unsigned n);
xt::xarray<double> eval_warp(int p, const xt::xarray<double> &xnodes,
                             const xt::xarray<double> &xout);
xt::xarray<double> eval_shift(int p, double pval, const xt::xarray<double> &l1,
                              const xt::xarray<double> &l2, const xt::xarray<double> &l3);
xt::xarray<double> generate_tetrahedron_nodes(unsigned n);
xt::xarray<double> warp_shift_face_3d(int p, double pval, double pval2,
                                      const xt::xarray<double> &l1, const xt::xarray<double> &l2,
                                      const xt::xarray<double> &l3, const xt::xarray<double> &l4);
xt::xarray<double> conversion_equilateral_xyz_to_rst(const xt::xarray<double> &coords);
xt::xarray<double> conversion_rst_to_abc(const xt::xarray<double> &rst);

xt::xarray<double> d_matrix_3d(const xt::xarray<double> &v, const xt::xarray<double> &gv);
xt::xarray<double> vandermonde_3d(unsigned n, const xt::xarray<double> &rst);
xt::xarray<double> vandermonde_3d_tensor(unsigned n, const xt::xarray<double> &rst);
xt::xarray<double> grad_simplexp_3d(const xt::xarray<double> &abc, int i, int j, int k);
xt::xarray<double> grad_tensorp_3d(const xt::xarray<double> &rst, int i, int j, int k);
xt::xarray<double> grad_vandermonde_3d(unsigned n, const xt::xarray<double> &rst);
xt::xarray<double> grad_vandermonde_3d_tensor(unsigned n, const xt::xarray<double> &rst);

xt::xarray<double> simplexp_3d(const xt::xarray<double> &abc, int i, int j, int k);
xt::xarray<double> tensorp_3d(const xt::xarray<double> &abc, int i, int j, int k);

}  // namespace oiseau::dg::nodal::utils
