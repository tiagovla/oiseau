// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace oiseau::mesh {

class Cell;
using CellType = const Cell *;

enum class CellKind {
  Undefined = 0,
  Point,
  Interval,
  Triangle,
  Quadrilateral,
  Tetrahedron,
  Hexahedron
};

CellType get_cell_type(const CellKind cell);

class Cell {
 protected:
  CellKind m_kind = CellKind::Undefined;
  std::string m_name;
  int m_dim{};
  std::vector<std::vector<std::vector<std::vector<int>>>> m_topology;
  std::pair<std::vector<double>, std::array<std::size_t, 2>> m_geometry;
  CellType m_facet = nullptr;
  CellType m_edge = nullptr;

 public:
  virtual ~Cell() = default;
  virtual CellType facet() { return m_facet; };
  virtual CellType edge() { return m_edge; };
  std::string_view name() const;
  int dimension() const;
  int num_entities(int dim) const;
  std::vector<std::vector<int>> get_entity_vertices(int dim) const;
  std::vector<std::vector<int>> get_sub_entities(int dim0, int dim1) const;
  const std::vector<std::vector<std::vector<std::vector<int>>>> &topology() const {
    return m_topology;
  }
  int num_sub_entities(int dim) const;
  CellKind kind() const { return m_kind; }
};

class PointCell : public Cell {
 public:
  PointCell();
};

class IntervalCell : public Cell {
 public:
  IntervalCell();
};

class TriangleCell : public Cell {
 private:
 public:
  TriangleCell();
};

class QuadrilateralCell : public Cell {
 public:
  QuadrilateralCell();
};

class TetrahedronCell : public Cell {
 public:
  TetrahedronCell();
};

class HexahedronCell : public Cell {
 public:
  HexahedronCell();
};

}  // namespace oiseau::mesh
