// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/io/gmsh.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oiseau/io/gmsh_file.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau::io {

namespace detail {
oiseau::mesh::CellType gmsh_celltype_to_oiseau_celltype(const std::size_t s) {
  static const std::unordered_map<std::size_t, oiseau::mesh::CellKind> gmsh_to_kind = {
      {15, oiseau::mesh::CellKind::Point},      {1, oiseau::mesh::CellKind::Interval},
      {2, oiseau::mesh::CellKind::Triangle},    {3, oiseau::mesh::CellKind::Quadrilateral},
      {4, oiseau::mesh::CellKind::Tetrahedron}, {5, oiseau::mesh::CellKind::Hexahedron},
  };
  auto it = gmsh_to_kind.find(s);
  if (it == gmsh_to_kind.end()) {
    throw std::runtime_error("Unknown Gmsh cell type: " + std::to_string(s));
  }

  return oiseau::mesh::get_cell_type(it->second);
}
}  // namespace detail

oiseau::mesh::Mesh gmsh_read_from_string(const std::string &content) {
  std::istringstream stream(content);
  return gmsh_read_from_stream(stream);
}

oiseau::mesh::Mesh gmsh_read_from_path(const std::filesystem::path &path) {
  std::ifstream f_handler(path);
  return gmsh_read_from_stream(f_handler);
}

oiseau::mesh::Mesh gmsh_read_from_stream(std::istream &f_handler) {
  GMSHFile file = GMSHFile(f_handler);
  std::vector<double> x;
  std::vector<std::vector<std::size_t>> conn;
  std::vector<oiseau::mesh::CellType> cell_types;

  x.reserve(file.nodes_section.num_nodes * 3);
  for (auto &block : file.nodes_section.blocks) {
    x.insert(x.end(), block.node_coords.begin(), block.node_coords.end());
  }

  cell_types.reserve(file.elements_section.num_elements);
  conn.reserve(file.elements_section.num_elements);

  for (const auto &block : file.elements_section.blocks) {
    std::size_t elem_size = block.data.size() / block.num_elements_in_block;
    for (std::size_t i = 0; i < block.num_elements_in_block; ++i) {
      std::vector<std::size_t> tmp;
      tmp.reserve(elem_size - 1);
      for (std::size_t j = 1; j < elem_size; ++j) {
        tmp.emplace_back(block.data[i * elem_size + j] - 1);
      }
      cell_types.emplace_back(detail::gmsh_celltype_to_oiseau_celltype(block.element_type));
      conn.emplace_back(std::move(tmp));
    }
  }

  oiseau::mesh::Geometry geometry = oiseau::mesh::Geometry(std::move(x), 3);
  oiseau::mesh::Topology topology = oiseau::mesh::Topology(std::move(conn), std::move(cell_types));
  oiseau::mesh::Mesh mesh(std::move(topology), std::move(geometry));
  return mesh;
};

void gmsh_write(const std::string &filename, const oiseau::mesh::Mesh &mesh) {
  throw std::logic_error("gmsh_write not implemented yet.");
};

}  // namespace oiseau::io
