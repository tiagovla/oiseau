// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// IWYU pragma: begin_exports
#include "oiseau/dg/connectivity.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/geometry.hpp"
#include "oiseau/dg/solver/linear_dg.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/time/integrator.hpp"
#include "oiseau/time/lsrk54.hpp"
#include "oiseau/time/rk4.hpp"
#include "oiseau/time/run.hpp"
// IWYU pragma: end_exports
