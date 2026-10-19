# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import meshio
import numpy as np

points = np.array(
    [
        # Hexahedron region (left)
        [0.0, 0.0, 0.0],  # 0
        [0.5, 0.0, 0.0],  # 1
        [0.5, 1.0, 0.0],  # 2
        [0.0, 1.0, 0.0],  # 3
        [0.0, 0.0, 1.0],  # 4
        [0.5, 0.0, 1.0],  # 5
        [0.5, 1.0, 1.0],  # 6
        [0.0, 1.0, 1.0],  # 7
        # Tetrahedron region (right)
        [1.0, 0.0, 0.0],  # 8
        [1.0, 1.0, 0.0],  # 9
        [1.0, 0.0, 1.0],  # 10
        [1.0, 1.0, 1.0],  # 11
    ]
)

# Cells
hexahedra = np.array([[0, 1, 2, 3, 4, 5, 6, 7]])
tetrahedra = np.array(
    [
        [1, 8, 2, 5],
        [5, 2, 8, 11],
    ]
)
cells = [
    ("hexahedron", hexahedra),
    ("tetra", tetrahedra),
]

# Tagging: hexa = geo 1, tetra = geo 2
cell_data = {
    "gmsh:physical": [
        np.array([1]),  # hexahedron in physical group 1
        np.array([2, 2]),  # tetrahedra in physical group 2
    ],
    "gmsh:geometrical": [
        np.array([1]),  # geometrical ID 1
        np.array([2, 2]),  # geometrical ID 2
    ],
}

# Point dim tags: (dim=3, tag=geo_id)
# Assign geo ID 1 to points 0-7, geo ID 2 to 8-11
point_data = {
    "gmsh:dim_tags": np.array(
        [
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 1),
            (3, 2),
            (3, 2),
            (3, 2),
            (3, 2),
        ],
        dtype=np.int32,
    )
}

mesh = meshio.Mesh(
    points=points, cells=cells, cell_data=cell_data, point_data=point_data
)

meshio.write(
    "mesh3d_tetra_hexa_block.msh",
    mesh,
    file_format="gmsh",
    binary=False,
)
