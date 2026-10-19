// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_triangle.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xoperation.hpp>
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

RefTriangle::RefTriangle(unsigned order) : RefElement(order) {
  this->m_np = ((order + 1) * (order + 2)) / 2;
  this->m_nfp = order + 1;
  this->m_r = detail::equilateral_xy_to_rs(detail::generate_triangle_nodes(this->m_order));
  this->m_v = this->vandermonde(this->m_r);
  this->m_gv = this->grad_vandermonde(this->m_r);
  this->m_d = this->grad_operator(this->m_v, this->m_gv);
}

xt::xarray<double> RefTriangle::basis_function(const xt::xarray<double> &ab, int i, int j) {
  xt::xarray<double> a = xt::col(ab, 0);
  xt::xarray<double> b = xt::col(ab, 1);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  auto h2 = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  return std::numbers::sqrt2 * h1 * h2 * xt::pow(1 - b, i);
}

xt::xarray<double> RefTriangle::grad_basis_function(const xt::xarray<double> &ab, int i, int j) {
  xt::xarray<double> a = xt::col(ab, 0);
  xt::xarray<double> b = xt::col(ab, 1);

  xt::xarray<double> fa = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  xt::xarray<double> dfa = oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, a);
  xt::xarray<double> gb = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  xt::xarray<double> dgb = oiseau::utils::grad_jacobi_p(j, 2.0 * i + 1.0, 0.0, b);

  xt::xarray<double> dmodedr = dfa * gb;
  xt::xarray<double> dmodeds = dfa * (gb * (0.5 * (1 + a)));

  xt::xarray<double> tmp = dgb * xt::pow(0.5 * (1 - b), i);

  if (i > 0) {
    dmodedr *= xt::pow(0.5 * (1 - b), i - 1);
    dmodeds *= xt::pow(0.5 * (1 - b), i - 1);
    tmp += -0.5 * i * gb * xt::pow(0.5 * (1 - b), i - 1);
  }

  dmodeds = dmodeds + fa * tmp;
  dmodedr = pow(2.0, i + 0.5) * dmodedr;
  dmodeds = pow(2.0, i + 0.5) * dmodeds;

  return xt::stack(xt::xtuple(dmodedr, dmodeds), 1);
}

xt::xarray<double> RefTriangle::vandermonde(const xt::xarray<double> &rs) const {
  const auto ab = detail::rs_to_ab(rs);

  const std::size_t n_points = rs.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)(order+2)/2

  xt::xarray<double> output = xt::zeros<double>({n_points, n_basis});

  std::size_t index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order - i; ++j, ++index) {
      xt::col(output, index) = this->basis_function(ab, i, j);
    }
  }
  return output;
}

xt::xarray<double> RefTriangle::grad_vandermonde(const xt::xarray<double> &rs) const {
  const auto ab = detail::rs_to_ab(rs);

  const std::size_t n_points = rs.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)(order+2)/2
  const std::size_t dimensions = 2;        // dr and ds

  xt::xarray<double> output = xt::zeros<double>({n_points, n_basis, dimensions});

  std::size_t index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order - i; ++j, ++index) {
      auto tmp = this->grad_basis_function(ab, i, j);
      xt::view(output, xt::all(), index, xt::all()) = tmp;
    }
  }
  return output;
}

xt::xarray<double> RefTriangle::grad_operator(const xt::xarray<double> &v,
                                              const xt::xarray<double> &gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);

  const auto vt = xt::transpose(v);
  const auto dvt_dr = xt::transpose(gvr);
  const auto dvt_ds = xt::transpose(gvs);

  const auto dr = xt::transpose(xt::linalg::solve(vt, dvt_dr));
  const auto ds = xt::transpose(xt::linalg::solve(vt, dvt_ds));

  return xt::stack(xt::xtuple(dr, ds), 2);
}

};  // namespace oiseau::dg::nodal

namespace oiseau::dg::nodal::detail {

xt::xarray<double> rs_to_ab(const xt::xarray<double> &rs) {
  auto r = xt::col(rs, 0);
  auto s = xt::col(rs, 1);
  auto ab = xt::zeros_like(rs);
  xt::xarray<double> denominator = 1.0 - s;
  xt::xarray<bool> is_one = xt::isclose(s, 1.0);  // handle s == 1 with tolerance
  xt::xarray<double> a = xt::where(is_one, -1.0, 2.0 * (1.0 + r) / denominator - 1.0);
  xt::col(ab, 0) = a;
  xt::col(ab, 1) = s;
  return ab;
}

xt::xarray<double> generate_triangle_nodes(unsigned order) {
  constexpr auto alp_opt =
      std::to_array({0.0000, 0.0000, std::numbers::sqrt2, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832,
                     1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258});
  double alpha = (order < 16) ? alp_opt[order - 1] : 5.0 / 3.0;
  std::size_t n_p = (order + 1) * (order + 2) / 2;
  xt::xarray<double> l1 = xt::zeros<double>({n_p});
  xt::xarray<double> l3 = xt::zeros<double>({n_p});
  double inv_n = 1.0 / order;

  std::size_t idx = 0;
  for (std::size_t i = 0; i < order + 1; ++i) {
    for (std::size_t j = 0; j < order + 1 - i; ++j, ++idx) {
      l1[idx] = i * inv_n;
      l3[idx] = j * inv_n;
    }
  }

  auto l2 = 1.0 - l1 - l3;
  auto blend1 = 4 * l2 * l3;
  auto blend2 = 4 * l1 * l3;
  auto blend3 = 4 * l1 * l2;

  auto e1 = -l2 + l3;
  auto e2 = -l3 + l1;
  auto e3 = -l1 + l2;

  auto warpf1 = detail::warp_factor(order, e1);
  auto warpf2 = detail::warp_factor(order, e2);
  auto warpf3 = detail::warp_factor(order, e3);

  auto warp1 = blend1 * warpf1 * xt::square(1 + alpha * l1);
  auto warp2 = blend2 * warpf2 * xt::square(1 + alpha * l2);
  auto warp3 = blend3 * warpf3 * xt::square(1 + alpha * l3);

  xt::xarray<double> out = xt::zeros<double>(xt::dynamic_shape<std::size_t>{n_p, 2});
  auto &&x = xt::col(out, 0);
  auto &&y = xt::col(out, 1);

  x = -l2 + l3;
  y = (-l2 - l3 + 2.0 * l1) / std::numbers::sqrt3;

  double pi = std::numbers::pi;
  x += warp1 + std::cos(2 * pi / 3) * warp2 + std::cos(4 * pi / 3) * warp3;
  y += std::sin(2 * pi / 3) * warp2 + std::sin(4 * pi / 3) * warp3;

  return out;
}

xt::xarray<double> equilateral_xy_to_rs(const xt::xarray<double> &xy) {
  constexpr double sqrt3 = std::numbers::sqrt3;
  auto x = xt::col(xy, 0);
  auto y = xt::col(xy, 1);
  auto l1 = (std::numbers::sqrt3 * y + 1.0) / 3.0;
  auto l2 = (-3.0 * x - sqrt3 * y + 2.0) / 6.0;
  auto l3 = (3.0 * x - sqrt3 * y + 2.0) / 6.0;

  xt::xarray<double> l_coords(xy.shape());
  xt::col(l_coords, 0) = -l2 + l3 - l1;
  xt::col(l_coords, 1) = -l2 - (l3 - l1);
  return l_coords;
}

xt::xarray<double> warp_factor(unsigned order, const xt::xarray<double> &rs) {
  auto lgl_r = oiseau::dg::nodal::utils::jacobi_gl(order, 0.0, 0.0);

  xt::xarray<double> r_eq = xt::linspace<double>(-1.0, 1.0, order + 1);
  xt::xarray<double> v_eq(xt::dynamic_shape<std::size_t>({r_eq.size(), order + 1}));
  for (size_t i = 0; i < order + 1; ++i) {
    xt::view(v_eq, xt::all(), i) = oiseau::utils::jacobi_p(i, 0.0, 0.0, r_eq);
  }

  xt::xarray<double> p_mat =
      xt::zeros<double>(xt::dynamic_shape<std::size_t>{order + 1, rs.size()});

  for (unsigned i = 0; i < order + 1; ++i) {
    xt::view(p_mat, i, xt::all()) = oiseau::utils::jacobi_p(i, 0.0, 0.0, rs);
  }

  auto l_mat = xt::linalg::solve(xt::transpose(v_eq), p_mat);
  auto warp = xt::linalg::dot(xt::transpose(l_mat), (-r_eq + lgl_r));
  auto zerof = xt::abs(rs) < 1.0 - 1e-10;
  auto sf = 1.0 - xt::pow((zerof * rs), 2);
  warp = warp / sf + warp * (zerof - 1);
  return warp;
}

}  // namespace oiseau::dg::nodal::detail
