# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import meshio


def main():
    # Define 8 cube vertices
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [1.0, 1.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [0.0, 0.0, 1.0],  # 4
            [1.0, 0.0, 1.0],  # 5
            [1.0, 1.0, 1.0],  # 6
            [0.0, 1.0, 1.0],  # 7
        ]
    )

    # Define 5 tetrahedra that fill the cube
    cells = [
        (
            "tetra",
            np.array(
                [
                    [0, 1, 3, 4],
                    [1, 2, 3, 6],
                    [1, 3, 4, 6],
                    [1, 4, 5, 6],
                    [3, 4, 6, 7],
                ]
            ),
        )
    ]

    # Create and write the mesh
    mesh = meshio.Mesh(points=points, cells=cells)
    meshio.write(
        "mesh3d_tetra_block.msh",
        mesh,
        file_format="gmsh",
        binary=False,
    )


if __name__ == "__main__":
    main()
