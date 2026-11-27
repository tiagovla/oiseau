// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/mesh.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/region.hpp"
#include "oiseau/mesh/topology.hpp"

using namespace oiseau::mesh;

Topology& Mesh::topology() { return _topology; }
const Topology& Mesh::topology() const { return _topology; }
Geometry& Mesh::geometry() { return _geometry; }
const Geometry& Mesh::geometry() const { return _geometry; }

void Mesh::add_region(std::string name, std::vector<std::size_t> cells) {
  if (name.empty()) throw std::invalid_argument("Region name cannot be empty");
  const auto cell_types = std::as_const(_topology).cell_types();
  if (cell_types.size() != _topology.n_cells()) {
    throw std::invalid_argument("Cell connectivity and cell type counts must match");
  }
  const auto dimension = _topology.dimension();
  for (const auto index : cells) {
    if (index >= cell_types.size()) {
      throw std::out_of_range("Region cell index is outside the mesh");
    }
    if (cell_types[index]->dimension() != dimension) {
      throw std::invalid_argument("Region cells must have the mesh dimension");
    }
  }

  const auto [region, inserted] =
      _regions.try_emplace(name, _regions.size() + 1, name, std::move(cells));
  if (!inserted) region->second.add_cells(std::move(cells));
}

const Region& Mesh::region(const std::string& name) const { return _regions.at(name); }
