// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <memory>
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
namespace oiseau::dg::nodal {

class Element {
 public:
  Element(std::shared_ptr<RefElement> ref_elem, xt::xarray<double> nodes);

  const RefElement& reference() const;
  unsigned order() const;
  const xt::xarray<double>& nodes() const;

 private:
  std::shared_ptr<RefElement> m_reference;
  xt::xarray<double> m_nodes;
};
}  // namespace oiseau::dg::nodal
