// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/cell.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <xtensor/containers/xadapt.hpp>

namespace oiseau::mesh {

CellType get_cell_type(const CellKind cell_kind) {
  static std::unordered_map<CellKind, std::unique_ptr<Cell>> cache;

  auto it = cache.find(cell_kind);
  if (it != cache.end()) {
    return it->second.get();
  }

  std::unique_ptr<Cell> cell;
  switch (cell_kind) {
  case CellKind::Point:
    cell = std::make_unique<PointCell>();
    break;
  case CellKind::Interval:
    cell = std::make_unique<IntervalCell>();
    break;
  case CellKind::Triangle:
    cell = std::make_unique<TriangleCell>();
    break;
  case CellKind::Quadrilateral:
    cell = std::make_unique<QuadrilateralCell>();
    break;
  case CellKind::Tetrahedron:
    cell = std::make_unique<TetrahedronCell>();
    break;
  case CellKind::Hexahedron:
    cell = std::make_unique<HexahedronCell>();
    break;
  default:
    throw std::runtime_error("Unknown cell type");
  }
  cache.emplace(cell_kind, std::move(cell));
  return cache[cell_kind].get();
}

std::string_view Cell::name() const { return m_name; }

int Cell::dimension() const { return m_dim; }

std::vector<std::vector<int>> Cell::get_sub_entities(int dim0, int dim1) const {
  return m_topology[dim0][dim1];
};

std::vector<std::vector<int>> Cell::get_entity_vertices(int dim) const {
  std::vector<std::vector<int>> slice;
  slice.reserve(m_topology[dim].size());
  std::transform(m_topology[dim].begin(), m_topology[dim].end(), std::back_inserter(slice),
                 [](const auto& s) { return s[0]; });
  return slice;
}

int Cell::num_sub_entities(int dim) const {
  if (dim <= m_dim) {
    return m_topology[dim].size();
  }
  return 0;
}

PointCell::PointCell() {
  m_name = "point";
  m_kind = CellKind::Point;
  m_dim = 0;
  // Reference-point vertex numbering:
  //
  //  0
  m_geometry = {
      {0.0},
      {1, 1},
  };
  m_topology = {
      {
          {{0}},
      },
  };
};

IntervalCell::IntervalCell() {
  m_name = "interval";
  m_kind = CellKind::Interval;
  m_dim = 1;
  // Reference-interval vertex and cell numbering:
  //
  //  0 -------- 1
  //       cell 0
  m_geometry = {
      {0.0, 1.0},
      {2, 1},
  };
  m_topology = {
      {
          {{0}, {0}},
          {{1}, {0}},
      },
      {
          {{0, 1}, {0}},
      },
  };
};

TriangleCell::TriangleCell() {
  m_name = "triangle";
  m_kind = CellKind::Triangle;
  m_dim = 2;
  // Reference-triangle vertex, edge, and face numbering:
  //
  //        2
  //       / \
  //    1 /  0\ 0
  //     /_____
  //    0   2   1
  //
  // Edges: 0=(1,2), 1=(0,2), 2=(0,1).
  // Face: 0=(0,1,2).
  m_geometry = {
      {0.0, 0.0, 1.0, 0.0, 0.0, 1.0},
      {3, 2},
  };
  m_topology = {
      {
          {{0}, {1, 2}, {0}},
          {{1}, {0, 2}, {0}},
          {{2}, {0, 1}, {0}},
      },
      {
          {{1, 2}, {0}, {0}},
          {{0, 2}, {1}, {0}},
          {{0, 1}, {2}, {0}},
      },
      {
          {{0, 1, 2}, {0, 1, 2}, {0}},
      },
  };
  m_facet = get_cell_type(CellKind::Interval);
  m_edge = get_cell_type(CellKind::Point);
}

QuadrilateralCell::QuadrilateralCell() {
  m_name = "quadrilateral";
  m_kind = CellKind::Quadrilateral;
  m_dim = 2;

  // Reference-quadrilateral vertex, edge, and face numbering:
  //
  //  3 ----2---- 2
  //  |           |
  //  3     0     1
  //  |           |
  //  0 ----0---- 1
  //
  // Edges: 0=(0,1), 1=(1,2), 2=(2,3), 3=(3,0).
  // Face: 0=(0,1,2,3).
  m_geometry = {
      {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0},
      {4, 2},
  };

  m_topology = {
      {
          {{0}, {0, 3}, {0}},
          {{1}, {0, 1}, {0}},
          {{2}, {1, 2}, {0}},
          {{3}, {2, 3}, {0}},
      },
      {
          {{0, 1}, {0, 1}, {0}},
          {{1, 2}, {1, 2}, {0}},
          {{2, 3}, {2, 3}, {0}},
          {{3, 0}, {3, 0}, {0}},
      },
      {
          {{0, 1, 2, 3}, {0, 1, 2, 3}, {0}},
      },
  };
  m_facet = get_cell_type(CellKind::Interval);
  m_edge = get_cell_type(CellKind::Point);
}

TetrahedronCell::TetrahedronCell() {
  m_name = "tetrahedron";
  m_kind = CellKind::Tetrahedron;
  m_dim = 3;
  // Reference-tetrahedron vertex, edge, and face numbering:
  //
  //          3
  //         /|\
  //        / | \
  //       1--|--2
  //        \ | /
  //         \|/
  //          0
  //
  // Edges: 0=(2,3), 1=(1,3), 2=(1,2), 3=(0,3), 4=(0,2), 5=(0,1).
  // Faces: 0=(1,2,3), 1=(0,2,3), 2=(0,1,3), 3=(0,1,2).
  m_geometry = {
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
      {4, 3},
  };
  m_topology = {
      {
          {{0}, {3, 4, 5}, {1, 2, 3}, {0}},
          {{1}, {1, 2, 5}, {0, 2, 3}, {0}},
          {{2}, {0, 2, 4}, {0, 1, 3}, {0}},
          {{3}, {0, 1, 3}, {0, 1, 2}, {0}},
      },
      {
          {{2, 3}, {0}, {0, 1}, {0}},
          {{1, 3}, {1}, {0, 2}, {0}},
          {{1, 2}, {2}, {0, 3}, {0}},
          {{0, 3}, {3}, {1, 2}, {0}},
          {{0, 2}, {4}, {1, 3}, {0}},
          {{0, 1}, {5}, {2, 3}, {0}},
      },
      {
          {{1, 2, 3}, {0, 1, 2}, {0}, {0}},
          {{0, 2, 3}, {0, 3, 4}, {1}, {0}},
          {{0, 1, 3}, {1, 3, 5}, {2}, {0}},
          {{0, 1, 2}, {2, 4, 5}, {3}, {0}},
      },
      {
          {{0, 1, 2, 3}, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3}, {0}},
      },
  };
  m_facet = get_cell_type(CellKind::Triangle);
  m_edge = get_cell_type(CellKind::Interval);
}

HexahedronCell::HexahedronCell() {
  m_name = "hexahedron";
  m_kind = CellKind::Hexahedron;
  m_dim = 3;

  // Reference-hexahedron vertex and edge numbering:
  //
  //        7 -------- 6
  //       /|         /|
  //      4 -------- 5 |
  //      | |        | |
  //      | 3 -------|-2
  //      |/         |/
  //      0 -------- 1
  //
  // Bottom edges: 0=(0,1), 1=(1,2), 2=(2,3), 3=(3,0).
  // Top edges:    4=(4,5), 5=(5,6), 6=(6,7), 7=(7,4).
  // Vertical:      8=(0,4), 9=(1,5), 10=(2,6), 11=(3,7).
  // Faces: 0=(0,1,2,3), 1=(4,5,6,7), 2=(0,1,5,4),
  //        3=(1,2,6,5), 4=(2,3,7,6), 5=(3,0,4,7).
  m_geometry = {
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0},
      {8, 3},
  };

  m_topology = {
      {
          {{0}, {0, 3, 8}, {0, 2, 5}, {0}},
          {{1}, {0, 1, 9}, {0, 2, 3}, {0}},
          {{2}, {1, 2, 10}, {0, 3, 4}, {0}},
          {{3}, {2, 3, 11}, {0, 4, 5}, {0}},
          {{4}, {4, 7, 8}, {1, 2, 5}, {0}},
          {{5}, {4, 5, 9}, {1, 2, 3}, {0}},
          {{6}, {5, 6, 10}, {1, 3, 4}, {0}},
          {{7}, {6, 7, 11}, {1, 4, 5}, {0}},
      },
      {
          {{0, 1}, {0}, {0, 2}, {0}},
          {{1, 2}, {1}, {0, 3}, {0}},
          {{2, 3}, {2}, {0, 4}, {0}},
          {{3, 0}, {3}, {0, 5}, {0}},
          {{4, 5}, {4}, {1, 2}, {0}},
          {{5, 6}, {5}, {1, 3}, {0}},
          {{6, 7}, {6}, {1, 4}, {0}},
          {{7, 4}, {7}, {1, 5}, {0}},
          {{0, 4}, {8}, {2, 5}, {0}},
          {{1, 5}, {9}, {2, 3}, {0}},
          {{2, 6}, {10}, {3, 4}, {0}},
          {{3, 7}, {11}, {4, 5}, {0}},
      },
      {
          {{0, 1, 2, 3}, {0, 1, 2, 3}, {0}, {0}},
          {{4, 5, 6, 7}, {4, 5, 6, 7}, {1}, {0}},
          {{0, 1, 5, 4}, {0, 4, 8, 9}, {2}, {0}},
          {{1, 2, 6, 5}, {1, 5, 9, 10}, {3}, {0}},
          {{2, 3, 7, 6}, {2, 6, 10, 11}, {4}, {0}},
          {{3, 0, 4, 7}, {3, 7, 8, 11}, {5}, {0}},
      },
      {
          {{0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7}, {0}},
      },
  };
  m_facet = get_cell_type(CellKind::Quadrilateral);
  m_edge = get_cell_type(CellKind::Interval);
}

}  // namespace oiseau::mesh
