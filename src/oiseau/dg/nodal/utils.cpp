// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/utils.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xoperation.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/utils/math.hpp"

namespace oiseau::dg::nodal::utils {

std::pair<xt::xarray<double>, xt::xarray<double>> jacobi_gq(unsigned n, double alpha, double beta) {
  if (n == 0) {
    xt::xarray<double> x = {(alpha - beta) / (alpha + beta + 2)};
    xt::xarray<double> w = {2};
    return make_pair(std::move(x), std::move(w));
  }
  xt::xarray<double> h1 = 2 * xt::linspace<double>(0, n, n + 1) + alpha + beta;
  auto h1s = xt::view(h1, xt::range(0, n));
  auto lins = xt::arange<double>(1, n + 1);
  xt::xarray<double> m =
      xt::diag(-1.0 / 2.0 * (std::pow(alpha, 2) - std::pow(beta, 2)) / (h1 + 2) / h1) +
      xt::diag(2. / (h1s + 2) *
                   xt::sqrt(lins * (lins + alpha + beta) * (lins + alpha) * (lins + beta) /
                            (h1s + 1) / (h1s + 3)),
               1);
  if (alpha + beta == 0) m(0, 0) = 0;
  m += xt::transpose(m);
  auto [d, v] = xt::linalg::eigh(m);
  double cte = std::pow(2, (alpha + beta + 1)) / (alpha + beta + 1) * std::tgamma(alpha + 1) *
               std::tgamma(beta + 1) / std::tgamma(alpha + beta + 1);
  xt::xarray<double> w = cte * xt::square(xt::view(v, 0, xt::all())), x = d;
  return {x, w};
}

xt::xarray<double> jacobi_gl(unsigned n, double alpha, double beta) {
  if (n == 1) {
    return xt::xarray<double>{-1.0, 1.0};
  }
  auto [x, _] = jacobi_gq(n - 2, alpha + 1.0, beta + 1.0);
  return xt::concatenate(xt::xtuple(-xt::ones<double>({1}), x, xt::ones<double>({1})));
}

xt::xarray<double> vandermonde_1d(unsigned n, const xt::xarray<double> &r) {
  std::array<size_t, 2> shape = {r.size(), n + 1};
  xt::xtensor<double, 2> output(shape);
  for (size_t i = 0; i < n + 1; ++i) {
    xt::view(output, xt::all(), i) = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);
  }
  return output;
}

xt::xarray<double> grad_vandermonde_1d(unsigned n, const xt::xarray<double> &r) {
  std::array<size_t, 2> shape = {r.size(), n + 1};
  auto output = xt::xtensor<double, 2>::from_shape(shape);
  for (size_t i = 0; i < n + 1; ++i) {
    xt::view(output, xt::all(), i) = oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, r);
  }
  return output;
}

xt::xarray<double> d_matrix_1d(const xt::xarray<double> &v, const xt::xarray<double> &gv) {
  return xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gv)));
}

xt::xarray<double> d_matrix_1d(unsigned n, const xt::xarray<double> &r) {
  auto v = vandermonde_1d(n, r);
  auto gv = grad_vandermonde_1d(n, r);
  return d_matrix_1d(v, gv);
}

xt::xarray<double> generate_triangle_equidistant_nodes(unsigned n) {
  std::size_t n_p = (n + 1) * (n + 2) / 2;
  xt::xarray<double> l1 = xt::zeros<double>({n_p});
  xt::xarray<double> l3 = xt::zeros<double>({n_p});
  xt::xarray<double> out = xt::zeros<double>(xt::xarray<double>::shape_type{n_p, 3});

  std::size_t idx = 0;
  for (std::size_t i = 0; i < n + 1; ++i) {
    for (std::size_t j = 0; j < n + 1 - i; ++j, ++idx) {
      l1[idx] = i * 1.0 / n;
      l3[idx] = j * 1.0 / n;
    }
  }

#define l2 (1.0 - l1 - l3)
  xt::view(out, xt::all(), 0) = -l2 + l3;
  xt::view(out, xt::all(), 1) = (-l2 - l3 + 2.0 * l1) / std::numbers::sqrt3;
#undef l2
  return out;
}

xt::xarray<double> warp_factor(unsigned n, const xt::xarray<double> &rout) {
  auto lgl_r = oiseau::dg::nodal::utils::jacobi_gl(n, 0.0, 0.0);
  auto r_eq = xt::linspace<double>(-1.0, 1.0, n + 1);
  auto v_eq = oiseau::dg::nodal::utils::vandermonde_1d(n, r_eq);
  xt::xarray<double> p_mat = xt::zeros<double>(xt::xarray<double>::shape_type{n + 1, rout.size()});
  for (unsigned i = 0; i < n + 1; ++i) {
    xt::view(p_mat, i, xt::all()) = oiseau::utils::jacobi_p(i, 0.0, 0.0, rout);
  }
  auto l_mat = xt::linalg::solve(xt::transpose(v_eq), p_mat);
  auto warp = xt::linalg::dot(xt::transpose(l_mat), (-r_eq + lgl_r));
  auto zerof = xt::abs(rout) < 1.0 - 1e-10;
  auto sf = 1.0 - xt::pow((zerof * rout), 2);
  warp = warp / sf + warp * (zerof - 1);
  return warp;
}

xt::xarray<double> conversion_equilateral_xy_to_rs(const xt::xarray<double> &coords) {
  auto l_coords = xt::zeros_like(coords);
  auto x = xt::col(coords, 0);
  auto y = xt::col(coords, 1);
  auto l1 = (std::numbers::sqrt3 * y + 1.0) / 3.0;
  auto l2 = (-3.0 * x - std::numbers::sqrt3 * y + 2.0) / 6.0;
  auto l3 = (3.0 * x - std::numbers::sqrt3 * y + 2.0) / 6.0;
  xt::col(l_coords, 0) = -l2 + l3 - l1;
  xt::col(l_coords, 1) = -l2 - (l3 - l1);
  return l_coords;
}

xt::xarray<double> simplexp_2d(const xt::xarray<double> &ab, int i, int j) {
  xt::xarray<double> a = xt::col(ab, 0);
  xt::xarray<double> b = xt::col(ab, 1);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  auto h2 = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  return std::numbers::sqrt2 * h1 * h2 * xt::pow(1 - b, i);
}

xt::xarray<double> tensorp_2d(const xt::xarray<double> &rs, int i, int j) {
  xt::xarray<double> r = xt::col(rs, 0);
  xt::xarray<double> s = xt::col(rs, 1);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);
  auto h2 = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);
  return h1 * h2;
}

xt::xarray<double> tensorp_3d(const xt::xarray<double> &rst, int i, int j, int k) {
  xt::xarray<double> r = xt::col(rst, 0);
  xt::xarray<double> s = xt::col(rst, 1);
  xt::xarray<double> t = xt::col(rst, 2);

  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);  // P_i(r)
  auto h2 = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);  // P_j(s)
  auto h3 = oiseau::utils::jacobi_p(k, 0.0, 0.0, t);  // P_k(t)

  return h1 * h2 * h3;
}

xt::xarray<double> grad_tensorp_2d(const xt::xarray<double> &rs, int i, int j) {
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

xt::xarray<double> grad_tensorp_3d(const xt::xarray<double> &rst, int i, int j, int k) {
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

xt::xarray<double> grad_simplexp_2d(const xt::xarray<double> &ab, int i, int j) {
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

xt::xarray<double> vandermonde_2d(unsigned n, const xt::xarray<double> &rs) {
  auto ab = conversion_rs_to_ab(rs);
  std::array<size_t, 2> shape = {rs.shape()[0], ((n + 1) * (n + 2)) / 2};
  xt::xtensor<double, 2> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n - i; ++j, ++index) {
      xt::col(output, index) = simplexp_2d(ab, i, j);
    }
  }
  return output;
}

xt::xarray<double> vandermonde_2d_tensor(unsigned n, const xt::xarray<double> &rs) {
  std::array<size_t, 2> shape = {rs.shape()[0], (n + 1) * (n + 1)};
  xt::xtensor<double, 2> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n; ++j, ++index) {
      xt::col(output, index) = xt::col(output, index) = tensorp_2d(rs, i, j);
    }
  }
  return output;
}

xt::xarray<double> grad_vandermonde_2d(unsigned n, const xt::xarray<double> &rs) {
  auto ab = conversion_rs_to_ab(rs);
  std::array<size_t, 3> shape = {rs.shape()[0], ((n + 1) * (n + 2)) / 2, 2};
  xt::xtensor<double, 3> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n - i; ++j, ++index) {
      auto tmp = grad_simplexp_2d(ab, i, j);
      xt::view(output, xt::all(), index, 0) = xt::col(tmp, 0);
      xt::view(output, xt::all(), index, 1) = xt::col(tmp, 1);
    }
  }
  return output;
}

xt::xarray<double> grad_vandermonde_2d_tensor(unsigned n, const xt::xarray<double> &rs) {
  std::array<size_t, 3> shape = {rs.shape()[0], (n + 1) * (n + 1), 2};
  xt::xtensor<double, 3> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n; ++j, ++index) {
      auto tmp = grad_tensorp_2d(rs, i, j);
      xt::view(output, xt::all(), index, 0) = xt::col(tmp, 0);
      xt::view(output, xt::all(), index, 1) = xt::col(tmp, 1);
    }
  }
  return output;
}

xt::xarray<double> d_matrix_2d(const xt::xarray<double> &v, const xt::xarray<double> &gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);
  auto dr = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvr)));
  auto ds = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvs)));
  return xt::stack(xt::xtuple(dr, ds), 2);
}

xt::xarray<double> conversion_rs_to_ab(const xt::xarray<double> &rs) {
  auto ab = xt::zeros_like(rs);
  auto r = xt::col(rs, 0);
  auto s = xt::col(rs, 1);
  for (unsigned i = 0; i < rs.shape()[0]; ++i) {
    ab(i, 0) = (s[i] == 1) ? -1 : 2 * (1 + r[i]) / (1 - s[i]) - 1;
  }
  xt::col(ab, 1) = s;
  return ab;
}

xt::xarray<double> generate_triangle_nodes(unsigned n) {
  constexpr auto alp_opt =
      std::to_array({0.0000, 0.0000, std::numbers::sqrt2, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832,
                     1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258});
  double alpha = (n < 16) ? alp_opt[n - 1] : 5.0 / 3.0;
  std::size_t n_p = (n + 1) * (n + 2) / 2;
  xt::xarray<double> l1 = xt::zeros<double>({n_p});
  xt::xarray<double> l3 = xt::zeros<double>({n_p});
  xt::xarray<double> out = xt::zeros<double>(xt::xarray<double>::shape_type{n_p, 2});
  std::size_t idx = 0;
  for (std::size_t i = 0; i < n + 1; ++i) {
    for (std::size_t j = 0; j < n + 1 - i; ++j, ++idx) {
      l1[idx] = i * 1.0 / n;
      l3[idx] = j * 1.0 / n;
    }
  }
  xt::xarray<double> l2 = (1.0 - l1 - l3);
  auto blend1 = 4 * l2 * l3;
  auto blend2 = 4 * l1 * l3;
  auto blend3 = 4 * l1 * l2;
  xt::xarray<double> e1 = -l2 + l3;
  xt::xarray<double> e2 = -l3 + l1;
  xt::xarray<double> e3 = -l1 + l2;
  auto warpf1 = oiseau::dg::nodal::utils::warp_factor(n, e1);
  auto warpf2 = oiseau::dg::nodal::utils::warp_factor(n, e2);
  auto warpf3 = oiseau::dg::nodal::utils::warp_factor(n, e3);
  auto warp1 = blend1 * warpf1 * xt::square(1 + (alpha * l1));
  auto warp2 = blend2 * warpf2 * xt::square(1 + (alpha * l2));
  auto warp3 = blend3 * warpf3 * xt::square(1 + (alpha * l3));
  auto &&x = xt::col(out, 0);
  auto &&y = xt::col(out, 1);
  auto pi = xt::numeric_constants<double>::PI;
  xt::col(out, 0) = -l2 + l3;
  xt::col(out, 1) = (-l2 - l3 + 2.0 * l1) / std::numbers::sqrt3;
  x = x + 1 * warp1 + std::cos(2 * pi / 3) * warp2 + std::cos(4 * pi / 3) * warp3;
  y = y + 0 * warp1 + std::sin(2 * pi / 3) * warp2 + std::sin(4 * pi / 3) * warp3;
  return out;
}

xt::xarray<double> generate_quadrilateral_nodes(unsigned order) {
  xt::xarray<double> r1d = jacobi_gl(order, 0.0, 0.0);
  unsigned Np1 = order + 1;
  unsigned Np_quad = Np1 * Np1;
  std::vector<std::size_t> shape = {static_cast<std::size_t>(Np_quad), 2ul};
  xt::xarray<double> rs_nodes = xt::zeros<double>(shape);

  unsigned k = 0;
  for (unsigned i = 0; i < Np1; ++i) {
    for (unsigned j = 0; j < Np1; ++j) {
      rs_nodes(k, 0) = r1d(j);
      rs_nodes(k, 1) = r1d(i);
      k++;
    }
  }
  return rs_nodes;
}

xt::xarray<double> generate_hexahedron_nodes(unsigned order) {
  xt::xarray<double> r1d = jacobi_gl(order, 0.0, 0.0);

  unsigned Np1 = order + 1;
  std::size_t Np_hex = Np1 * Np1 * Np1;

  std::vector<std::size_t> shape = {Np_hex, 3ul};
  xt::xarray<double> rst_nodes = xt::zeros<double>(shape);

  unsigned k_idx = 0;
  for (unsigned k_loop = 0; k_loop < Np1; ++k_loop) {
    for (unsigned i_loop = 0; i_loop < Np1; ++i_loop) {
      for (unsigned j_loop = 0; j_loop < Np1; ++j_loop) {
        rst_nodes(k_idx, 0) = r1d(j_loop);
        rst_nodes(k_idx, 1) = r1d(i_loop);
        rst_nodes(k_idx, 2) = r1d(k_loop);
        k_idx++;
      }
    }
  }
  return rst_nodes;
}

xt::xarray<double> generate_tetrahedron_equidistant_nodes(unsigned n) {
  auto n_p = ((n + 1) * (n + 2) * (n + 3)) / 6;
  auto shape = xt::xarray<double>::shape_type{n_p, 3};
  xt::xarray<double> output = xt::zeros<double>(shape);
  std::size_t sk = 0;
  for (std::size_t i = 0; i < n + 1; ++i) {
    for (std::size_t j = 0; j < n + 1 - i; ++j) {
      for (std::size_t k = 0; k < n + 1 - i - j; ++k, ++sk) {
        output(sk, 0) = -1 + k * 2.0 / n;
        output(sk, 1) = -1 + j * 2.0 / n;
        output(sk, 2) = -1 + i * 2.0 / n;
      }
    }
  }
  return output;
}

xt::xarray<double> eval_warp(int p, const xt::xarray<double> &xnodes,
                             const xt::xarray<double> &xout) {
  xt::xarray<double> warp = xt::zeros_like(xout);
  xt::xarray<double> xeq = xt::zeros<double>({p + 1});
  for (int i = 0; i <= p; i++) xeq(i) = -1.0 + 2.0 * (p - i) / p;
  for (int i = 0; i <= p; i++) {
    xt::xarray<double> d = xt::ones_like(xout) * (xnodes(i) - xeq(i));
    for (int j = 1; j < p; j++) {
      if (i != j) d = d * (xout - xeq(j)) / (xeq(i) - xeq(j));
    }
    if (i != 0) d = -d / (xeq(i) - xeq(0));
    if (i != p) d = d / (xeq(i) - xeq(p));
    warp += d;
  }
  return warp;
}

xt::xarray<double> eval_shift(int p, double pval, const xt::xarray<double> &l1,
                              const xt::xarray<double> &l2, const xt::xarray<double> &l3) {
  auto gaussX = -jacobi_gl(p, 0.0, 0.0);
  auto blend1 = l2 * l3;
  auto blend2 = l1 * l3;
  auto blend3 = l1 * l2;
  auto warpfactor1 = 4 * eval_warp(p, gaussX, l3 - 1 * l2);
  auto warpfactor2 = 4 * eval_warp(p, gaussX, l1 - 1 * l3);
  auto warpfactor3 = 4 * eval_warp(p, gaussX, l2 - 1 * l1);
  auto warp1 = blend1 * warpfactor1 * (1 + xt::pow(pval * l1, 2));
  auto warp2 = blend2 * warpfactor2 * (1 + xt::pow(pval * l2, 2));
  auto warp3 = blend3 * warpfactor3 * (1 + xt::pow(pval * l3, 2));
  auto pi = xt::numeric_constants<double>::PI;
  auto dx = 1 * warp1 + cos(2.0 * pi / 3.0) * warp2 + cos(4.0 * pi / 3.0) * warp3;
  auto dy = 0 * warp1 + sin(2.0 * pi / 3.0) * warp2 + sin(4.0 * pi / 3.0) * warp3;
  return xt::stack(xt::xtuple(dx, dy), 1);
}

xt::xarray<double> warp_shift_face_3d(int p, double pval, const xt::xarray<double> &l2,
                                      const xt::xarray<double> &l3, const xt::xarray<double> &l4) {
  return eval_shift(p, pval, l2, l3, l4);
}

xt::xarray<double> generate_tetrahedron_nodes(unsigned p) {
  xt::xarray<double> alphastore = {
      0,      0,       0,      0.1002, 1.1332, 1.5608, 1.3413, 1.2577,
      1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655,
  };
  double alpha = p > 14 ? 1 : alphastore[p];
  double tol = 1e-10;
  auto rst = generate_tetrahedron_equidistant_nodes(p);
  auto r = xt::col(rst, 0), s = xt::col(rst, 1), t = xt::col(rst, 2);
  xt::xarray<double> l1 = (1 + t) / 2, l2 = (1 + s) / 2, l3 = -(1 + r + s + t) / 2,
                     l4 = (1 + r) / 2;
  xt::xtensor<double, 1> v1 = {-1.0, -1.0 / std::numbers::sqrt3, -1.0 / std::sqrt(6)};
  xt::xtensor<double, 1> v2 = {1.0, -1.0 / std::numbers::sqrt3, -1.0 / std::sqrt(6)};
  xt::xtensor<double, 1> v3 = {0.0, 2.0 / std::numbers::sqrt3, -1.0 / std::sqrt(6)};
  xt::xtensor<double, 1> v4 = {0.0, 0.0, 3.0 / std::sqrt(6)};
  xt::xtensor<double, 2> t1 = xt::zeros<double>({4, 3});
  xt::xtensor<double, 2> t2 = xt::zeros<double>({4, 3});
  xt::row(t1, 0) = v2 - 1 * v1;
  xt::row(t1, 1) = v2 - 1 * v1;
  xt::row(t1, 2) = v3 - 1 * v2;
  xt::row(t1, 3) = v3 - 1 * v1;
  xt::row(t2, 0) = v3 - 0.5 * (v1 + v2);
  xt::row(t2, 1) = v4 - 0.5 * (v1 + v2);
  xt::row(t2, 2) = v4 - 0.5 * (v2 + v3);
  xt::row(t2, 3) = v4 - 0.5 * (v1 + v3);
  for (int i = 0; i < 4; ++i) {
    xt::row(t1, i) /= xt::linalg::norm(xt::row(t1, i));
    xt::row(t2, i) /= xt::linalg::norm(xt::row(t2, i));
  }
  auto xyz = xt::linalg::outer(l3, v1) + xt::linalg::outer(l4, v2) + xt::linalg::outer(l2, v3) +
             xt::linalg::outer(l1, v4);
  xt::xarray<double> shift = xt::zeros<double>(xyz.shape());
  for (int face = 0; face < 4; ++face) {
    xt::xarray<double> la, lb, lc, ld;
    if (face == 0) {
      la = l1;
      lb = l2;
      lc = l3;
      ld = l4;
    } else if (face == 1) {
      la = l2;
      lb = l1;
      lc = l3;
      ld = l4;
    } else if (face == 2) {
      la = l3;
      lb = l1;
      lc = l4;
      ld = l2;
    } else if (face == 3) {
      la = l4;
      lb = l1;
      lc = l3;
      ld = l2;
    }
    auto warp = warp_shift_face_3d(p, alpha, lb, lc, ld);
    auto warp1 = xt::col(warp, 0);
    auto warp2 = xt::col(warp, 1);
    xt::xarray<double> blend = lb * lc * ld;
    xt::xarray<double> denom = (lb + 0.5 * la) * (lc + 0.5 * la) * (ld + 0.5 * la);
    auto ids = xt::flatten_indices(xt::argwhere(denom > tol));
    xt::view(blend, xt::keep(ids)) = (1 + xt::pow(alpha * xt::view(la, xt::keep(ids)), 2)) *
                                     xt::view(blend, xt::keep(ids)) /
                                     xt::view(denom, xt::keep(ids));
    shift = shift + xt::linalg::outer(blend * warp1, xt::row(t1, face)) +
            xt::linalg::outer(blend * warp2, xt::row(t2, face));
    auto mask = (la < tol) && (((lb > tol) + (lc > tol) + (ld > tol)) < 3);
    ids = xt::flatten_indices(xt::argwhere(mask));
    xt::view(shift, xt::keep(ids), xt::all()) =
        xt::linalg::outer(xt::view(warp1, xt::keep(ids)), xt::row(t1, face)) +
        xt::linalg::outer(xt::view(warp2, xt::keep(ids)), xt::row(t2, face));
  }
  xyz = xyz + shift;
  return xyz;
}

xt::xarray<double> conversion_equilateral_xyz_to_rst(const xt::xarray<double> &coords) {
  const double sqrt3 = std::numbers::sqrt3;
  const double sqrt6 = std::sqrt(6.0);
  xt::xtensor<double, 1> v1 = {-1, -1 / sqrt3, -1 / sqrt6};
  xt::xtensor<double, 1> v2 = {1, -1 / sqrt3, -1 / sqrt6};
  xt::xtensor<double, 1> v3 = {0, 2 / sqrt3, -1 / sqrt6};
  xt::xtensor<double, 1> v4 = {0, 0 / sqrt3, 3 / sqrt6};
  auto rhs = xt::transpose(coords - 0.5 * (v2 + v3 + v4 - v1));
  auto lhs = 0.5 * xt::stack(xt::xtuple(v2 - v1, v3 - v1, v4 - v1), 1);
  auto rst = xt::linalg::solve(lhs, rhs);
  return xt::transpose(rst);
}

xt::xarray<double> simplexp_3d(const xt::xarray<double> &abc, int i, int j, int k) {
  xt::xarray<double> a = xt::col(abc, 0);
  xt::xarray<double> b = xt::col(abc, 1);
  xt::xarray<double> c = xt::col(abc, 2);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  auto h2 = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  auto h3 = oiseau::utils::jacobi_p(k, 2.0 * (i + j) + 2.0, 0.0, c);
  return 2 * std::numbers::sqrt2 * h1 * h2 * xt::pow(1 - b, i) * h3 * xt::pow(1 - c, i + j);
}

xt::xarray<double> conversion_rst_to_abc(const xt::xarray<double> &rst) {
  auto abc = xt::zeros_like(rst);
  auto &&r = xt::col(rst, 0);
  auto &&s = xt::col(rst, 1);
  auto &&t = xt::col(rst, 2);
  auto &&a = xt::col(abc, 0);
  auto &&b = xt::col(abc, 1);
  auto &&c = xt::col(abc, 2);
  for (unsigned n = 0; n < rst.shape()[0]; ++n) {
    if ((s[n] + t[n]) != 0) {
      a[n] = 2 * ((1 + r[n]) / (-s[n] - t[n])) - 1;
    } else {
      a[n] = -1;
    }
    if (t[n] != 1) {
      b[n] = 2 * (1 + s[n]) / (1 - t[n]) - 1;
    } else {
      b[n] = -1;
    }
  }
  c = t;
  return abc;
}

xt::xarray<double> vandermonde_3d(unsigned n, const xt::xarray<double> &rst) {
  auto abc = conversion_rst_to_abc(rst);
  std::array<size_t, 2> shape = {rst.shape()[0], ((n + 1) * (n + 2) * (n + 3)) / 6};
  xt::xtensor<double, 2> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n - i; ++j) {
      for (unsigned k = 0; k <= n - i - j; ++k, ++index) {
        xt::col(output, index) = simplexp_3d(abc, i, j, k);
      }
    }
  }
  return output;
}

xt::xarray<double> vandermonde_3d_tensor(unsigned n, const xt::xarray<double> &rst) {
  std::array<size_t, 2> shape = {rst.shape()[0], (n + 1) * (n + 1) * (n + 1)};
  xt::xtensor<double, 2> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n; ++j) {
      for (unsigned k = 0; k <= n; ++k, ++index) {
        xt::col(output, index) = tensorp_3d(rst, i, j, k);
      }
    }
  }
  return output;
}

xt::xarray<double> grad_simplexp_3d(const xt::xarray<double> &abc, int i, int j, int k) {
  xt::xarray<double> a = xt::col(abc, 0);
  xt::xarray<double> b = xt::col(abc, 1);
  xt::xarray<double> c = xt::col(abc, 2);
  xt::xarray<double> fa = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  xt::xarray<double> dfa = oiseau::utils::grad_jacobi_p(i, 0.0, 0.0, a);
  xt::xarray<double> gb = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  xt::xarray<double> dgb = oiseau::utils::grad_jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  xt::xarray<double> hc = oiseau::utils::jacobi_p(k, 2.0 * (i + j) + 2.0, 0.0, c);
  xt::xarray<double> dhc = oiseau::utils::grad_jacobi_p(k, 2.0 * (i + j) + 2.0, 0.0, c);

  xt::xarray<double> v3dr, v3ds, v3dt, tmp;

  v3dr = dfa * (gb * hc);
  if (i > 0) v3dr = v3dr * xt::pow(0.5 * (1 - b), i - 1);
  if (i + j > 0) v3dr = v3dr * xt::pow(0.5 * (1 - c), i + j - 1);

  v3ds = 0.5 * (1 + a) * v3dr;
  tmp = dgb * xt::pow(0.5 * (1 - b), i);
  if (i > 0) tmp = tmp + (-0.5 * i) * gb * xt::pow(0.5 * (1 - b), i - 1);
  if (i + j > 0) tmp = tmp * xt::pow(0.5 * (1 - c), i + j - 1);
  tmp = fa * (tmp * hc);
  v3ds = v3ds + tmp;

  v3dt = 0.5 * (1 + a) * v3dr + 0.5 * (1 + b) * tmp;
  tmp = dhc * xt::pow(0.5 * (1 - c), i + j);
  if (i + j > 0) tmp = tmp - 0.5 * (i + j) * hc * xt::pow(0.5 * (1 - c), (i + j - 1));
  tmp = fa * gb * tmp;
  tmp = tmp * xt::pow(0.5 * (1 - b), i);
  v3dt = v3dt + tmp;

  v3dr = v3dr * std::pow(2, 2 * i + j + 1.5);
  v3ds = v3ds * std::pow(2, 2 * i + j + 1.5);
  v3dt = v3dt * std::pow(2, 2 * i + j + 1.5);

  return xt::stack(xt::xtuple(v3dr, v3ds, v3dt), 1);
}

xt::xarray<double> grad_vandermonde_3d(unsigned n, const xt::xarray<double> &rst) {
  auto abc = conversion_rst_to_abc(rst);
  std::array<size_t, 3> shape = {rst.shape()[0], ((n + 1) * (n + 2) * (n + 3)) / 6, 3};
  xt::xtensor<double, 3> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n - i; ++j) {
      for (unsigned k = 0; k <= n - i - j; ++k, ++index) {
        auto tmp = grad_simplexp_3d(abc, i, j, k);
        xt::view(output, xt::all(), index, 0) = xt::col(tmp, 0);
        xt::view(output, xt::all(), index, 1) = xt::col(tmp, 1);
        xt::view(output, xt::all(), index, 2) = xt::col(tmp, 2);
      }
    }
  }
  return output;
}

xt::xarray<double> grad_vandermonde_3d_tensor(unsigned n, const xt::xarray<double> &rst) {
  std::array<size_t, 3> shape = {rst.shape()[0], (n + 1) * (n + 1) * (n + 1), 3};
  xt::xtensor<double, 3> output(shape);
  int index = 0;
  for (unsigned i = 0; i <= n; ++i) {
    for (unsigned j = 0; j <= n; ++j) {
      for (unsigned k = 0; k <= n; ++k, ++index) {
        auto tmp = grad_tensorp_3d(rst, i, j, k);
        xt::view(output, xt::all(), index, 0) = xt::col(tmp, 0);
        xt::view(output, xt::all(), index, 1) = xt::col(tmp, 1);
        xt::view(output, xt::all(), index, 2) = xt::col(tmp, 2);
      }
    }
  }
  return output;
}

xt::xarray<double> d_matrix_3d(const xt::xarray<double> &v, const xt::xarray<double> &gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);
  auto gvt = xt::view(gv, xt::all(), xt::all(), 2);
  auto dr = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvr)));
  auto ds = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvs)));
  auto dt = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvt)));
  return xt::stack(xt::xtuple(dr, ds, dt), 2);
}

}  // namespace oiseau::dg::nodal::utils
