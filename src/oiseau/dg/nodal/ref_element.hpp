// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <memory>
#include <stdexcept>
#include <xtensor/containers/xarray.hpp>

#include "xtensor/core/xtensor_forward.hpp"

namespace oiseau::dg::nodal {

enum class RefElementType { Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron };

class RefElement {
 public:
  virtual ~RefElement() = default;

  inline const xt::xarray<double>& v() const { return m_v; }
  inline const xt::xarray<double>& gv() const { return m_gv; }
  inline const xt::xarray<double>& d() const { return m_d; }
  inline const xt::xarray<double>& r() const { return m_r; }

  inline unsigned order() const { return m_order; }
  inline unsigned number_of_nodes() const { return m_np; }
  inline unsigned number_of_face_nodes() const { return m_nfp; }

  virtual xt::xarray<double> vandermonde(const xt::xarray<double>&) const = 0;

 protected:
  explicit RefElement(unsigned order) : m_order(order) {
    if (order == 0) throw std::invalid_argument("Order must be greater than 0");
  }

  unsigned m_order;
  unsigned m_np{};
  unsigned m_nfp{};
  xt::xarray<double> m_v;
  xt::xarray<double> m_gv;
  xt::xarray<double> m_d;
  xt::xarray<double> m_r;
};

std::shared_ptr<RefElement> get_ref_element(RefElementType type, unsigned order);

}  // namespace oiseau::dg::nodal
