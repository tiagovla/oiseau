# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import oiseau.mesh

triangle = oiseau.mesh.get_cell_type("triangle")
vertices = triangle.get_entity_vertices(1)
print(vertices)
