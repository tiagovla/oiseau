// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/io/gmsh_file.hpp"

#include <array>
#include <cstddef>
#include <istream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <xtensor/containers/xadapt.hpp>

enum { PREFIX = '$' };

namespace oiseau::io {

namespace detail {
std::size_t gmsh_nodes_per_cell(const std::size_t s) {
  switch (s) {
  case 1:
    return 2;  // line
  case 2:
    return 3;  // triangle
  case 3:
    return 4;  // quad
  case 4:
    return 4;  // tetra
  case 5:
    return 8;  // hexahedron
  case 6:
    return 6;  // prism?
  case 7:
    return 5;  // wedge
  case 15:
    return 1;  // vertex
  default:
    throw std::runtime_error("Unknown GMSH cell type: " + std::to_string(s));
  }
}

template <class T>
concept Readable = std::is_convertible_v<T, std::string_view> || std::is_arithmetic_v<T>;

template <Readable T, int N>
std::array<T, N> from_file(std::istream& f, bool is_binary = false) {
  std::array<T, N> arr{};
  if (is_binary) {
    f.read(reinterpret_cast<char*>(arr.data()), sizeof(T) * N);
  } else {
    for (auto& value : arr) f >> value;
  }
  return arr;
}

template <Readable T>
std::vector<T> from_file(std::istream& f, std::size_t n, bool is_binary = false) {
  std::vector<T> vec(n);
  if (is_binary) {
    f.read(reinterpret_cast<char*>(vec.data()), sizeof(T) * n);
  } else {
    for (std::size_t i = 0; i < n; i++) f >> vec.at(i);
  }
  return vec;
}

MeshFormatSection mesh_format_handler(std::istream& f_handler) {
  auto [version] = from_file<double, 1>(f_handler);
  if (version != 4.1) {
    throw std::runtime_error(
        "Unsupported GMSH version detected."
        "Please ensure you are using version 4.1.");
  }
  auto [is_binary] = from_file<int, 1>(f_handler);
  auto [size_t_size] = from_file<std::size_t, 1>(f_handler);
  if (is_binary) {
    f_handler.get();  // skip NF
    auto [verify_one] = from_file<int, 1>(f_handler, true);
    if (verify_one != 1) throw std::runtime_error("Invalid GMSH file");
    if (size_t_size != sizeof(std::size_t)) throw std::runtime_error("Invalid GMSH file");
  }
  return {version, is_binary, size_t_size};
}

PhysicalNamesSection physical_names_handler(std::istream& f_handler) {
  auto [num_phys_names] = from_file<int, 1>(f_handler);
  std::vector<int> dimensions;
  std::vector<int> physical_tags;
  std::vector<std::string> names;
  dimensions.reserve(num_phys_names);
  physical_tags.reserve(num_phys_names);
  names.reserve(num_phys_names);
  for (int i = 0; i < num_phys_names; i++) {
    auto [dim, tag] = from_file<int, 2>(f_handler);
    auto [name] = from_file<std::string, 1>(f_handler);
    dimensions.emplace_back(dim);
    physical_tags.emplace_back(tag);
    names.emplace_back(name);
  }
  return {num_phys_names, std::move(dimensions), std::move(physical_tags), std::move(names)};
}

EntitiesSection entities_handler(std::istream& f_handler, bool is_binary) {
  auto quantity = from_file<std::size_t, 4>(f_handler, is_binary);
  std::array<std::vector<EntityEntry>, 4> blocks;
  for (std::size_t i = 0; i < 4; i++) blocks.at(i).reserve(quantity.at(i));
  for (int d = 0; d < 4; d++) {
    for (std::size_t j = 0; j < quantity[d]; j++) {
      auto [tag] = from_file<int, 1>(f_handler, is_binary);
      auto bounding_coods = from_file<double>(f_handler, (d == 0) ? 3 : 6, is_binary);
      auto [num_physicals] = from_file<std::size_t, 1>(f_handler, is_binary);
      auto physical_tags = from_file<int>(f_handler, num_physicals, is_binary);
      if (d > 0) {
        auto [num_BREP] = from_file<std::size_t, 1>(f_handler, is_binary);
        auto bounding_entities = from_file<int>(f_handler, num_BREP, is_binary);
        blocks.at(d).emplace_back(tag, std::move(bounding_coods), std::move(physical_tags),
                                  std::move(bounding_entities));
      } else {
        blocks.at(d).emplace_back(tag, std::move(bounding_coods), std::move(physical_tags));
      }
    }
  }
  return EntitiesSection(std::move(blocks));
}

NodesSection nodes_handler(std::istream& f_handler, bool is_binary) {
  auto [num_entity_blocks, total_num_nodes, min_node_tag, max_node_tag] =
      from_file<std::size_t, 4>(f_handler, is_binary);

  std::vector<NodesBlock> blocks;
  blocks.reserve(num_entity_blocks);

  for (std::size_t i = 0; i < num_entity_blocks; i++) {
    auto [dim, entity_tag, parametric] = from_file<int, 3>(f_handler, is_binary);
    auto [quantity] = from_file<std::size_t, 1>(f_handler, is_binary);
    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords;
    node_tags.reserve(quantity);
    node_coords.reserve(quantity * 3);

    for (std::size_t j = 0; j < quantity; j++) {
      auto [tag] = from_file<std::size_t, 1>(f_handler, is_binary);
      node_tags.emplace_back(tag);
    }
    for (std::size_t j = 0; j < quantity; j++) {
      auto xyz = from_file<double, 3>(f_handler, is_binary);
      node_coords.insert(node_coords.end(), xyz.begin(), xyz.end());
    }
    blocks.emplace_back(dim, entity_tag, parametric, quantity, std::move(node_tags),
                        std::move(node_coords));
  }
  return {num_entity_blocks, total_num_nodes, min_node_tag, max_node_tag, std::move(blocks)};
}

ElementSection elements_handler(std::istream& f_handler, bool is_binary) {
  auto [num_element_blocks, num_elements, min_element_tag, max_element_tag] =
      from_file<std::size_t, 4>(f_handler, is_binary);
  std::vector<ElementBlock> blocks;
  blocks.reserve(num_element_blocks);
  for (std::size_t i = 0; i < num_element_blocks; i++) {
    auto [entity_dim, entity_tag, element_type] = from_file<int, 3>(f_handler, is_binary);
    auto [num_elements_in_block] = from_file<std::size_t, 1>(f_handler, is_binary);
    auto size = gmsh_nodes_per_cell(element_type);
    auto tmp_conn =
        from_file<std::size_t>(f_handler, (1 + size) * num_elements_in_block, is_binary);
    blocks.emplace_back(entity_dim, entity_tag, element_type, num_elements_in_block,
                        std::move(tmp_conn));
  }
  return {num_element_blocks, num_elements, min_element_tag, max_element_tag, std::move(blocks)};
};

void skip_to_end_of_environment(std::istream& f_handler) {
  std::string line;
  while (getline(f_handler, line)) {
    if (line.starts_with(PREFIX)) break;
  }
}
}  // namespace detail

void GMSHFile::read(std::istream& f_handler) {
  using namespace detail;
  if (f_handler.fail()) throw std::runtime_error("Could not read file stream");
  std::string line;
  bool is_binary;
  while (getline(f_handler, line)) {
    if (line.starts_with(PREFIX)) {
      line = line.substr(1);
      if (line == "MeshFormat") {
        mesh_format_section = mesh_format_handler(f_handler);
        is_binary = mesh_format_section.is_binary;
      } else if (line == "PhysicalNames") {
        physical_names_section = physical_names_handler(f_handler);
      } else if (line == "Entities") {
        entities_section = entities_handler(f_handler, is_binary);
      } else if (line == "Nodes") {
        nodes_section = nodes_handler(f_handler, is_binary);
      } else if (line == "Elements") {
        elements_section = elements_handler(f_handler, is_binary);
      }
      skip_to_end_of_environment(f_handler);
    }
  }
}
}  // namespace oiseau::io
