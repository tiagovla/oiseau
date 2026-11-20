# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import oiseau.cpp as cpp

CellKind = cpp.mesh.CellKind
Mesh = cpp.mesh.Mesh


def get_cell_type(cell_type: str | CellKind):
    if isinstance(cell_type, str):
        cell_type = CellKind.__members__[cell_type.capitalize()]
    return cpp.mesh.get_cell_type(cell_type)
