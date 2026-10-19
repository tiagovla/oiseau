// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <istream>
#include <string>
#include <utility>
#include <vector>

namespace oiseau::io {

struct MeshFormatSection {
  double version{};
  int is_binary{};
  std::size_t data_size{};

  MeshFormatSection() = default;
  MeshFormatSection(double version, int is_binary, std::size_t data_size)
      : version(version), is_binary(is_binary), data_size(data_size) {}
};

struct EntityEntry {
  std::size_t tag;
  std::vector<double> boundary_coords;
  std::vector<int> physical_tags;
  std::vector<int> bounding_tags;

  EntityEntry(std::size_t tag, std::vector<double>&& boundary_coords,
              std::vector<int>&& physical_tags, std::vector<int>&& bounding_tags)
      : tag(tag),
        boundary_coords(std::move(boundary_coords)),
        physical_tags(std::move(physical_tags)),
        bounding_tags(std::move(bounding_tags)) {}
  EntityEntry(std::size_t tag, std::vector<double>&& boundary_coords,
              std::vector<int>&& physical_tags)
      : tag(tag),
        boundary_coords(std::move(boundary_coords)),
        physical_tags(std::move(physical_tags)) {}
};

struct EntitiesSection {
  std::array<std::vector<EntityEntry>, 4> blocks{};

  EntitiesSection() = default;
  explicit EntitiesSection(std::array<std::vector<EntityEntry>, 4>&& blocks) : blocks(blocks) {}
};

struct PhysicalNamesSection {
  int num_physical_names{};
  std::vector<int> dimensions{};
  std::vector<int> physical_tags{};
  std::vector<std::string> names{};

  PhysicalNamesSection() = default;
  PhysicalNamesSection(int num_physical_names, std::vector<int>&& dimensions,
                       std::vector<int>&& physical_tags, std::vector<std::string>&& names)
      : num_physical_names(num_physical_names),
        dimensions(std::move(dimensions)),
        physical_tags(std::move(physical_tags)),
        names(std::move(names)) {};
};

struct NodesBlock {
  int entity_dim;
  int entity_tag;
  int parametric;
  std::size_t num_nodes_in_block;
  std::vector<std::size_t> node_tags;
  std::vector<double> node_coords;

  NodesBlock(int entity_dim, int entity_tag, int parametric, std::size_t num_nodes_in_block,
             std::vector<std::size_t>&& node_tags, std::vector<double>&& node_coords)
      : entity_dim(entity_dim),
        entity_tag(entity_tag),
        parametric(parametric),
        num_nodes_in_block(num_nodes_in_block),
        node_tags(std::move(node_tags)),
        node_coords(std::move(node_coords)) {}
};

struct NodesSection {
  std::size_t num_entity_blocks{};
  std::size_t num_nodes{};
  std::size_t min_node_tag{};
  std::size_t max_node_tag{};
  std::vector<NodesBlock> blocks{};

  NodesSection() = default;
  NodesSection(std::size_t num_entity_blocks, std::size_t num_nodes, std::size_t min_node_tag,
               std::size_t max_node_tag, std::vector<NodesBlock>&& blocks)
      : num_entity_blocks(num_entity_blocks),
        num_nodes(num_nodes),
        min_node_tag(min_node_tag),
        max_node_tag(max_node_tag),
        blocks(std::move(blocks)) {}
};

struct ElementBlock {
  int entity_dim;
  int entity_tag;
  int element_type;
  std::size_t num_elements_in_block;
  std::vector<std::size_t> data;

  ElementBlock(int entity_dim, int entity_tag, int element_type, std::size_t num_elements_in_block,
               std::vector<std::size_t>&& data)
      : entity_dim(entity_dim),
        entity_tag(entity_tag),
        element_type(element_type),
        num_elements_in_block(num_elements_in_block),
        data(std::move(data)) {}
};

struct ElementSection {
  std::size_t num_element_blocks{};
  std::size_t num_elements{};
  std::size_t min_element_tag{};
  std::size_t max_element_tag{};
  std::vector<ElementBlock> blocks{};

  ElementSection() = default;
  ElementSection(std::size_t num_element_blocks, std::size_t num_elements,
                 std::size_t min_element_tag, std::size_t max_element_tag,
                 std::vector<ElementBlock>&& blocks)
      : num_element_blocks(num_element_blocks),
        num_elements(num_elements),
        min_element_tag(min_element_tag),
        max_element_tag(max_element_tag),
        blocks(std::move(blocks)) {}
};

class GMSHFile {
 public:
  GMSHFile(MeshFormatSection&& mesh_format_section, PhysicalNamesSection&& physical_names_section,
           EntitiesSection&& entities_section, NodesSection&& nodes_section,
           ElementSection&& elements_section)
      : mesh_format_section(std::move(mesh_format_section)),
        physical_names_section(std::move(physical_names_section)),
        entities_section(std::move(entities_section)),
        nodes_section(std::move(nodes_section)),
        elements_section(std::move(elements_section)) {}

  explicit GMSHFile(std::istream& f_handler) { read(f_handler); }

  MeshFormatSection mesh_format_section;
  PhysicalNamesSection physical_names_section;
  EntitiesSection entities_section;
  NodesSection nodes_section;
  ElementSection elements_section;

 private:
  void read(std::istream& f_handler);
};

namespace detail {
MeshFormatSection mesh_format_handler(std::istream& f_handler);
PhysicalNamesSection physical_names_handler(std::istream& f_handler);
EntitiesSection entities_handler(std::istream& f_handler, bool is_binary);
NodesSection nodes_handler(std::istream& f_handler, bool is_binary);
ElementSection elements_handler(std::istream& f_handler, bool is_binary);

}  // namespace detail
}  // namespace oiseau::io
