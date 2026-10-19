// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <cstddef>
#include <span>
#include <vector>

#include "oiseau/mesh/cell.hpp"

namespace oiseau::mesh {

class Topology {
 public:
  Topology();
  Topology(std::vector<std::vector<std::size_t>>&& conn, std::vector<CellType>&& cell_types);
  Topology(Topology&&) = default;
  Topology(const Topology&) = default;
  Topology& operator=(Topology&&) = default;
  Topology& operator=(const Topology&) = default;
  ~Topology();
  std::span<CellType> cell_types();
  std::span<const CellType> cell_types() const;
  std::span<std::vector<std::size_t>> conn();
  std::span<const std::vector<std::size_t>> conn() const;
  std::span<std::vector<std::size_t>> e_to_e();
  std::span<std::vector<std::size_t>> e_to_f();
  std::size_t n_cells() const;
  int dimension() const;
  void calculate_connectivity();

 private:
  std::vector<std::vector<std::size_t>> m_conn;
  std::vector<std::vector<std::size_t>> m_e_to_v;
  std::vector<std::vector<std::size_t>> m_e_to_e;
  std::vector<std::vector<std::size_t>> m_e_to_f;
  std::vector<CellType> m_cell_types;
  mutable int m_dimension = -1;
};

}  // namespace oiseau::mesh
