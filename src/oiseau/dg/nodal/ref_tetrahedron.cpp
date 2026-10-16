// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_tetrahedron.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
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

RefTetrahedron::RefTetrahedron(unsigned order) : RefElement(order) {
  this->m_np = ((order + 1) * (order + 2) * (order + 3)) / 6;
  this->m_nfp = ((order + 1) * (order + 2)) / 2;
  this->m_r = detail::equilateral_xyz_to_rst(detail::generate_tetrahedron_nodes(this->m_order));
  this->m_v = RefTetrahedron::vandermonde(this->m_r);
  this->m_gv = this->grad_vandermonde(this->m_r);
  this->m_d = this->grad_operator(this->m_v, this->m_gv);
}

xt::xarray<double> RefTetrahedron::basis_function(const xt::xarray<double> &abc, int i, int j,
                                                  int k) const {
  xt::xarray<double> a = xt::col(abc, 0);
  xt::xarray<double> b = xt::col(abc, 1);
  xt::xarray<double> c = xt::col(abc, 2);
  auto h1 = oiseau::utils::jacobi_p(i, 0.0, 0.0, a);
  auto h2 = oiseau::utils::jacobi_p(j, 2.0 * i + 1.0, 0.0, b);
  auto h3 = oiseau::utils::jacobi_p(k, 2.0 * (i + j) + 2.0, 0.0, c);
  return 2 * std::numbers::sqrt2 * h1 * h2 * xt::pow(1 - b, i) * h3 * xt::pow(1 - c, i + j);
}

xt::xarray<double> RefTetrahedron::grad_basis_function(const xt::xarray<double> &abc, int i, int j,
                                                       int k) const {
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

xt::xarray<double> RefTetrahedron::vandermonde(const xt::xarray<double> &rst) const {
  auto abc = detail::rst_to_abc(rst);
  const std::size_t n_points = rst.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)(order+2)*(order+3)/6
  //
  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis}));

  std::size_t index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order - i; ++j) {
      for (unsigned k = 0; k <= this->m_order - i - j; ++k, ++index) {
        xt::col(output, index) = this->basis_function(abc, i, j, k);
      }
    }
  }
  return output;
}

xt::xarray<double> RefTetrahedron::grad_vandermonde(const xt::xarray<double> &rst) const {
  auto abc = detail::rst_to_abc(rst);
  const std::size_t n_points = rst.shape()[0];
  const std::size_t n_basis = this->m_np;  // (order+1)(order+2)*(order+3)/6
  //
  xt::xarray<double> output(xt::dynamic_shape<std::size_t>({n_points, n_basis, 3}));

  std::size_t index = 0;
  for (unsigned i = 0; i <= this->m_order; ++i) {
    for (unsigned j = 0; j <= this->m_order - i; ++j) {
      for (unsigned k = 0; k <= this->m_order - i - j; ++k, ++index) {
        auto tmp = this->grad_basis_function(abc, i, j, k);
        xt::view(output, xt::all(), index, xt::all()) = tmp;
      }
    }
  }
  return output;
}
//
xt::xarray<double> RefTetrahedron::grad_operator(const xt::xarray<double> &v,
                                                 const xt::xarray<double> &gv) const {
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

xt::xarray<double> rst_to_abc(const xt::xarray<double> &rst) {
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
  auto gaussX = -oiseau::dg::nodal::utils::jacobi_gl(p, 0.0, 0.0);
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

xt::xarray<double> equilateral_xyz_to_rst(const xt::xarray<double> &coords) {
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

}  // namespace oiseau::dg::nodal::detail
