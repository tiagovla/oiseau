// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_quadrilateral.hpp"

#include <cstddef>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xshape.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/utils.hpp"
#include "oiseau/utils/math.hpp"

namespace oiseau::dg::nodal {

RefQuadrilateral::RefQuadrilateral(unsigned order) : RefElement(order) {
  this->m_np = (order + 1) * (order + 1);
  this->m_nfp = order + 1;
  this->m_r = detail::generate_quadrilateral_nodes(this->m_order);
  this->m_v = RefQuadrilateral::vandermonde(this->m_r);
  this->m_gv = this->grad_vandermonde(this->m_r);
  this->m_d = grad_operator(this->m_v, this->m_gv);
}

xt::xarray<double> RefQuadrilateral::basis_function(const xt::xarray<double> &rs, int i, int j) {
  xt::xarray<double> r = xt::col(rs, 0);
  xt::xarray<double> s = xt::col(rs, 1);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);
  auto h2 = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);
  return h1 * h2;
}
//
xt::xarray<double> RefQuadrilateral::grad_basis_function(const xt::xarray<double> &rs, int i,
                                                         int j) {
  xt::xarray<double> r = xt::col(rs, 0);
  xt::xarray<double> s = xt::col(rs, 1);

  auto pi_r = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);        // P_i(r)
  auto dpi_r = oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, r);  // dP_i/dr
  auto pj_s = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);        // P_j(s)
  auto dpj_s = oiseau::utils::grad_jacobi_p(j, 0.0, 0.0, s);  // dP_j/ds

  auto dphidr = dpi_r * pj_s;  // ∂φ/∂r
  auto dphids = pi_r * dpj_s;  // ∂φ/∂s

  return xt::stack(xt::xtuple(dphidr, dphids), 1);  // shape: [N, 2]
}

xt::xarray<double> RefQuadrilateral::vandermonde(const xt::xarray<double> &rs) const {
  const std::size_t n_points = rs.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)*(order+1)

  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis}));
  int index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order; ++j, ++index) {
      xt::col(output, index) = this->basis_function(rs, i, j);
    }
  }
  return output;
}

xt::xarray<double> RefQuadrilateral::grad_vandermonde(const xt::xarray<double> &rs) const {
  const std::size_t n_points = rs.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)*(order+1)

  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis, 2}));

  int index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order; ++j, ++index) {
      auto tmp = this->grad_basis_function(rs, i, j);
      xt::view(output, xt::all(), index, xt::all()) = tmp;
    }
  }
  return output;
}

xt::xarray<double> RefQuadrilateral::grad_operator(const xt::xarray<double> &v,
                                                   const xt::xarray<double> &gv) const {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);

  const auto vt = xt::transpose(v);
  const auto dvt_dr = xt::transpose(gvr);
  const auto dvt_ds = xt::transpose(gvs);

  const auto dr = xt::transpose(xt::linalg::solve(vt, dvt_dr));
  const auto ds = xt::transpose(xt::linalg::solve(vt, dvt_ds));

  return xt::stack(xt::xtuple(dr, ds), 2);
}

}  // namespace oiseau::dg::nodal

namespace oiseau::dg::nodal::detail {

xt::xarray<double> generate_quadrilateral_nodes(unsigned order) {
  xt::xarray<double> r1d = utils::jacobi_gl(order, 0.0, 0.0);
  auto [R, S] = xt::meshgrid(r1d, r1d);
  return xt::stack(xt::xtuple(xt::flatten(R), xt::flatten(S)), 1);
}
//

}  // namespace oiseau::dg::nodal::detail
