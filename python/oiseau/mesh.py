# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import oiseau.cpp as cpp


def get_cell_type(cell_type: str):
    return cpp.mesh.get_cell_type(cell_type)
