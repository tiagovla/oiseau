// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <utility>

#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau::mesh {
class Mesh {
 public:
  Mesh() : _topology(), _geometry() {}
  Mesh(Topology &topology, Geometry &geometry) : _topology(topology), _geometry(geometry) {}
  Mesh(Topology &&topology, Geometry &&geometry)
      : _topology(std::move(topology)), _geometry(std::move(geometry)) {}

  Mesh(Mesh &&) = default;
  Mesh(const Mesh &) = default;
  Mesh &operator=(Mesh &&) = default;
  Mesh &operator=(const Mesh &) = default;
  ~Mesh() = default;

  Topology &topology();
  const Topology &topology() const;
  Geometry &geometry();
  const Geometry &geometry() const;

 private:
  Topology _topology;
  Geometry _geometry;
};

}  // namespace oiseau::mesh
