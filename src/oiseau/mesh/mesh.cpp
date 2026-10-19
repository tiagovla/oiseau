// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/mesh/mesh.hpp"

#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/topology.hpp"

using namespace oiseau::mesh;

Topology& Mesh::topology() { return _topology; }
const Topology& Mesh::topology() const { return _topology; }
Geometry& Mesh::geometry() { return _geometry; }
const Geometry& Mesh::geometry() const { return _geometry; }
