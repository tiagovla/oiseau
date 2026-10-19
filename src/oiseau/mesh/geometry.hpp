// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <cstddef>
#include <span>
#include <vector>

namespace oiseau::mesh {
class Geometry {
 public:
  Geometry();
  Geometry(std::vector<double>&& x, unsigned dim);
  Geometry(Geometry&&) = default;
  Geometry(const Geometry&) = default;
  Geometry& operator=(Geometry&&) = default;
  Geometry& operator=(const Geometry&) = default;
  ~Geometry();

  std::span<double> x();
  std::span<const double> x() const;
  std::span<double> x_at(std::size_t pos);
  unsigned dim() const;

 private:
  std::vector<double> m_x;
  unsigned m_dim = 3;
};
}  // namespace oiseau::mesh
