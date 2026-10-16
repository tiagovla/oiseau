// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_line.hpp"

#include <cstddef>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xshape.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/utils/integration.hpp"
#include "oiseau/utils/math.hpp"

namespace oiseau::dg::nodal {

RefLine::RefLine(unsigned order) : RefElement(order) {
  this->m_np = order + 1;
  this->m_nfp = 1;
  this->m_r = detail::generate_line_nodes(this->m_order);
  this->m_v = RefLine::vandermonde(this->m_r);
  this->m_gv = this->grad_vandermonde(this->m_r);
  this->m_d = this->grad_operator(this->m_v, this->m_gv);
}

xt::xarray<double> RefLine::basis_function(const xt::xarray<double>& r, int i) {
  return oiseau::utils::jacobi_p(i, 0.0, 0.0, r);
}

xt::xarray<double> RefLine::grad_basis_function(const xt::xarray<double>& r, int i) {
  return oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, r);
}

xt::xarray<double> RefLine::vandermonde(const xt::xarray<double>& r) const {
  const auto n_points = r.shape()[0];
  const auto n_basis = this->m_np;  // order + 1
  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis}));
  for (std::size_t i = 0; i <= this->m_order; ++i) {
    xt::view(output, xt::all(), i) = this->basis_function(r, i);
  }
  return output;
}

xt::xarray<double> RefLine::grad_vandermonde(const xt::xarray<double>& r) const {
  const auto n_points = r.shape()[0];
  const auto n_basis = this->m_np;  // order + 1
  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis}));
  for (size_t i = 0; i <= this->m_order; ++i) {
    xt::view(output, xt::all(), i) = grad_basis_function(r, i);
  }
  return output;
}

xt::xarray<double> RefLine::grad_operator(const xt::xarray<double>& v,
                                          const xt::xarray<double>& gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);

  const auto vt = xt::transpose(v);
  const auto dvt_dr = xt::transpose(gvr);
  const auto dr = xt::transpose(xt::linalg::solve(vt, dvt_dr));
  return dr;
}

}  // namespace oiseau::dg::nodal

namespace oiseau::dg::nodal::detail {

xt::xarray<double> generate_line_nodes(unsigned order) {
  if (order == 1) {
    return xt::xarray<double>{-1.0, 1.0};
  }
  auto [x, _] = oiseau::utils::integration::jacobi_gq(order - 2, 1.0, 1.0);
  return xt::concatenate(xt::xtuple(-xt::ones<double>({1}), x, xt::ones<double>({1})));
}
}  // namespace oiseau::dg::nodal::detail
