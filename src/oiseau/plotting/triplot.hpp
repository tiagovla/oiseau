// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "oiseau/mesh/mesh.hpp"
#include "oiseau/plotting/pyplot.hpp"

namespace oiseau::plotting {
void triplot(plt::AxesSubPlot &ax, oiseau::mesh::Mesh &mesh);
}
