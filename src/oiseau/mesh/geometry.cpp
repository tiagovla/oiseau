// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/geometry.hpp"

#include <cstddef>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace oiseau::mesh;

Geometry::Geometry() = default;
Geometry::~Geometry() = default;
std::span<double> Geometry::x() { return m_x; };
std::span<const double> Geometry::x() const { return m_x; };
std::span<double> Geometry::x_at(std::size_t pos) { return {&m_x[pos * m_dim], m_dim}; };
Geometry::Geometry(std::vector<double>&& x, unsigned dim) : m_x(std::move(x)), m_dim(dim) {};

Geometry::Geometry(std::initializer_list<std::initializer_list<double>> points) {
  const auto dimension = points.size() == 0 ? 0 : points.begin()->size();
  if (dimension == 0 || !std::in_range<unsigned>(dimension)) {
    throw std::invalid_argument("Geometry requires a positive, representable coordinate dimension");
  }
  m_dim = static_cast<unsigned>(dimension);
  for (const auto& point : points) {
    if (point.size() != dimension) {
      throw std::invalid_argument("Geometry points must all have the same coordinate dimension");
    }
    m_x.insert_range(m_x.end(), point);
  }
}

unsigned Geometry::dim() const { return m_dim; };
