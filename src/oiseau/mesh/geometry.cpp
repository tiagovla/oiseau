// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/geometry.hpp"

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

using namespace oiseau::mesh;

Geometry::Geometry() = default;
Geometry::~Geometry() = default;
std::span<double> Geometry::x() { return m_x; };
std::span<const double> Geometry::x() const { return m_x; };
std::span<double> Geometry::x_at(std::size_t pos) { return {&m_x[pos * m_dim], m_dim}; };
Geometry::Geometry(std::vector<double> &&x, unsigned dim) : m_x(std::move(x)), m_dim(dim) {};
unsigned Geometry::dim() const { return m_dim; };
