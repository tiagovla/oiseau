// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/ref_element.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/misc/xmanipulation.hpp>

#include "oiseau/dg/nodal/ref_hexahedron.hpp"
#include "oiseau/dg/nodal/ref_line.hpp"
#include "oiseau/dg/nodal/ref_quadrilateral.hpp"
#include "oiseau/dg/nodal/ref_tetrahedron.hpp"
#include "oiseau/dg/nodal/ref_triangle.hpp"

namespace oiseau::dg::nodal {

xt::xarray<double> RefElement::mass() const {
  const auto inverse = inv_v();
  return xt::linalg::dot(xt::transpose(inverse), inverse);
}

xt::xarray<double> RefElement::inv_mass() const { return xt::linalg::dot(m_v, xt::transpose(m_v)); }

xt::xarray<double> RefElement::inv_v() const { return xt::linalg::inv(m_v); }

std::shared_ptr<RefElement> get_ref_element(RefElementType type, unsigned order) {
  using Key = std::pair<RefElementType, unsigned>;

  struct KeyHash {
    std::size_t operator()(const Key& k) const {
      return std::hash<int>()(static_cast<int>(k.first)) ^ std::hash<unsigned>()(k.second);
    }
  };

  static std::unordered_map<Key, std::shared_ptr<RefElement>, KeyHash> cache;

  Key key{type, order};
  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }

  std::shared_ptr<RefElement> elem;
  switch (type) {
  case RefElementType::Line:
    elem = std::make_shared<RefLine>(order);
    break;
  case RefElementType::Triangle:
    elem = std::make_shared<RefTriangle>(order);
    break;
  case RefElementType::Quadrilateral:
    elem = std::make_shared<RefQuadrilateral>(order);
    break;
  case RefElementType::Tetrahedron:
    elem = std::make_shared<RefTetrahedron>(order);
    break;
  case RefElementType::Hexahedron:
    elem = std::make_shared<RefHexahedron>(order);
    break;
  default:
    throw std::invalid_argument("Unknown element type");
  }

  cache[key] = elem;
  return elem;
}

}  // namespace oiseau::dg::nodal
