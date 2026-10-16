// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_hexahedron.hpp"

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

RefHexahedron::RefHexahedron(unsigned order) : RefElement(order) {
  this->m_np = (order + 1) * (order + 1) * (order + 1);
  this->m_nfp = (order + 1) * (order + 1);
  this->m_r = detail::generate_hexahedron_nodes(this->m_order);
  this->m_v = RefHexahedron::vandermonde(this->m_r);
  this->m_gv = this->grad_vandermonde(this->m_r);
  this->m_d = this->grad_operator(this->m_v, this->m_gv);
}

xt::xarray<double> RefHexahedron::basis_function(const xt::xarray<double> &rst, int i, int j,
                                                 int k) {
  xt::xarray<double> r = xt::col(rst, 0);
  xt::xarray<double> s = xt::col(rst, 1);
  xt::xarray<double> t = xt::col(rst, 2);

  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);  // P_i(r)
  auto h2 = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);  // P_j(s)
  auto h3 = oiseau::utils::jacobi_p(k, 0.0, 0.0, t);  // P_k(t)

  return h1 * h2 * h3;
}

xt::xarray<double> RefHexahedron::grad_basis_function(const xt::xarray<double> &rst, int i, int j,
                                                      int k) {
  xt::xarray<double> r = xt::col(rst, 0);
  xt::xarray<double> s = xt::col(rst, 1);
  xt::xarray<double> t = xt::col(rst, 2);

  auto Pi_r = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);  // P_i(r)
  auto Pj_s = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);  // P_j(s)
  auto Pk_t = oiseau::utils::jacobi_p(k, 0.0, 0.0, t);  // P_k(t)

  auto dPi_r = oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, r);  // dP_i/dr
  auto dPj_s = oiseau::utils::grad_jacobi_p(j, 0.0, 0.0, s);  // dP_j/ds
  auto dPk_t = oiseau::utils::grad_jacobi_p(k, 0.0, 0.0, t);  // dP_k/dt

  auto dphidr = dPi_r * Pj_s * Pk_t;
  auto dphids = Pi_r * dPj_s * Pk_t;
  auto dphidt = Pi_r * Pj_s * dPk_t;

  return xt::stack(xt::xtuple(dphidr, dphids, dphidt), 1);  // shape: [N, 3]
}

xt::xarray<double> RefHexahedron::vandermonde(const xt::xarray<double> &rst) const {
  const std::size_t n_points = rst.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)*(order+1)*(order+1)

  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis}));
  int index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order; ++j) {
      for (unsigned k = 0; k <= this->m_order; ++k, ++index) {
        xt::col(output, index) = this->basis_function(rst, i, j, k);
      }
    }
  }
  return output;
}

xt::xarray<double> RefHexahedron::grad_vandermonde(const xt::xarray<double> &rst) const {
  const std::size_t n_points = rst.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)*(order+1)*(order+1)

  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis, 3}));
  std::size_t index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order; ++j) {
      for (unsigned k = 0; k <= this->m_order; ++k, ++index) {
        auto tmp = this->grad_basis_function(rst, i, j, k);
        xt::view(output, xt::all(), index, xt::all()) = tmp;
      }
    }
  }
  return output;
}

xt::xarray<double> RefHexahedron::grad_operator(const xt::xarray<double> &v,
                                                const xt::xarray<double> &gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);
  auto gvt = xt::view(gv, xt::all(), xt::all(), 2);

  const auto vt = xt::transpose(v);
  const auto dvt_dr = xt::transpose(gvr);
  const auto dvt_ds = xt::transpose(gvs);
  const auto dvt_dt = xt::transpose(gvt);

  const auto dr = xt::transpose(xt::linalg::solve(vt, dvt_dr));
  const auto ds = xt::transpose(xt::linalg::solve(vt, dvt_ds));
  const auto dt = xt::transpose(xt::linalg::solve(vt, dvt_dt));

  return xt::stack(xt::xtuple(dr, ds, dt), 2);
}

}  // namespace oiseau::dg::nodal

namespace oiseau::dg::nodal::detail {
xt::xarray<double> generate_hexahedron_nodes(unsigned order) {
  xt::xarray<double> r1d = utils::jacobi_gl(order, 0.0, 0.0);
  auto [R, S, T] = xt::meshgrid(r1d, r1d, r1d);
  return xt::stack(xt::xtuple(xt::flatten(R), xt::flatten(S), xt::flatten(T)), 1);
}

}  // namespace oiseau::dg::nodal::detail
