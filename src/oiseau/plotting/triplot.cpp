// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/plotting/triplot.hpp"

#include <cstddef>
#include <limits>
#include <vector>

#include "oiseau/mesh/mesh.hpp"

namespace oiseau::plotting {

void triplot(plt::AxesSubPlot& ax, oiseau::mesh::Mesh& mesh) {
  const auto& topology = mesh.topology();
  const auto& geometry = mesh.geometry();
  const auto connectivity = topology.conn();
  const auto coordinates = geometry.x();
  const auto dimension = geometry.dim();
  const auto cell_types = topology.cell_types();
  std::vector<double> x, y;
  constexpr auto separator = std::numeric_limits<double>::quiet_NaN();

  for (std::size_t i = 0; i < connectivity.size(); ++i) {
    const auto& vertices = connectivity[i];
    for (const auto& face_vertices : cell_types[i]->get_entity_vertices(1)) {
      for (const auto local_vertex : face_vertices) {
        const auto vertex = vertices[local_vertex];
        x.push_back(coordinates[vertex * dimension]);
        y.push_back(coordinates[vertex * dimension + 1]);
      }
      x.push_back(separator);
      y.push_back(separator);
    }
  }

  ax.plot(x, y, "color"_a = "k");
}

}  // namespace oiseau::plotting
