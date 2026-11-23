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

xt::xarray<double> tensorp_3d(const xt::xarray<double>& rst, int i, int j, int k) {
  xt::xarray<double> r = xt::col(rst, 0);
  xt::xarray<double> s = xt::col(rst, 1);
  xt::xarray<double> t = xt::col(rst, 2);

  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, r);  // P_i(r)
  auto h2 = oiseau::utils::jacobi_p(j, 0.0, 0.0, s);  // P_j(s)
  auto h3 = oiseau::utils::jacobi_p(k, 0.0, 0.0, t);  // P_k(t)

  return h1 * h2 * h3;
}

xt::xarray<double> grad_tensorp_3d(const xt::xarray<double>& rst, int i, int j, int k) {
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

xt::xarray<double> vandermonde_3d_tensor(unsigned n, const xt::xarray<double>& rst) {
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

xt::xarray<double> grad_vandermonde_3d_tensor(unsigned n, const xt::xarray<double>& rst) {
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

xt::xarray<double> d_matrix_3d(const xt::xarray<double>& v, const xt::xarray<double>& gv) {
  auto gvr = xt::view(gv, xt::all(), xt::all(), 0);
  auto gvs = xt::view(gv, xt::all(), xt::all(), 1);
  auto gvt = xt::view(gv, xt::all(), xt::all(), 2);
  auto dr = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvr)));
  auto ds = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvs)));
  auto dt = xt::transpose(xt::linalg::solve(xt::transpose(v), xt::transpose(gvt)));
  return xt::stack(xt::xtuple(dr, ds, dt), 2);
}

}  // namespace oiseau::dg::nodal::utils
