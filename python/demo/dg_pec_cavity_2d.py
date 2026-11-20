# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

from argparse import ArgumentParser
from math import pi, sin, sqrt

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize

from oiseau.dg import DGSpace, EquationSystem, Field, LinearDG2D, dt, dx, dy, run_rk4
from oiseau.mesh import CellKind, Mesh


def pec_wall(interior, _normal):
    ez, hx, hy = interior
    return [-ez, hx, hy]


def plot_ez(group_nodes, ez_groups):
    amplitude = max(np.abs(values).max() for values in ez_groups)
    norm = Normalize(vmin=-amplitude, vmax=amplitude)
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for nodes, values in zip(group_nodes, ez_groups, strict=True):
        # Triangulate each cell separately to preserve jumps across DG faces.
        for cell_nodes, ez_values in zip(nodes, values, strict=True):
            x, y = cell_nodes.T
            field_plot = ax.tripcolor(x, y, ez_values, shading="gouraud", cmap="RdBu_r", norm=norm)

    ax.set(aspect="equal", xlabel="x", ylabel="y", title=r"$E_z$ after half a period")
    ax.margins(0)
    fig.colorbar(field_plot, ax=ax, shrink=0.8, label=r"$E_z$")
    return fig


def main(output=None):
    mesh = Mesh(
        [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)],
        [[0, 1, 4, 3], [1, 2, 4], [2, 5, 4]],
        [CellKind.Quadrilateral, CellKind.Triangle, CellKind.Triangle],
    )
    space = DGSpace(mesh, order=[4, 6, 6])

    ez, hx, hy = Field("Ez"), Field("Hx"), Field("Hy")
    equations = EquationSystem(
        [ez, hx, hy],
        [dt(ez) == dx(hy) - dy(hx), dt(hx) == -dy(ez), dt(hy) == dx(ez)],
    )
    cavity = LinearDG2D(mesh, space, equations, wave_speed=1.0, boundary_state=pec_wall)
    initial = cavity.initial_state(ez, lambda x, y: sin(pi * x / 2) * sin(pi * y))

    half_period = 2 / sqrt(5)
    steps = 256
    state = run_rk4(cavity, initial, 0, half_period / steps, steps)

    # The exact electric field changes sign after half a period.
    field = equations.field_index(ez)
    ez_groups = [values[..., field] for values in state.groups]
    max_error = np.max(
        [
            np.abs(values + start[..., field]).max()
            for values, start in zip(ez_groups, initial.groups, strict=True)
        ]
    )
    print(f"Maximum Ez error after half a period: {max_error:.6f}")
    if not np.isfinite(max_error):
        raise SystemExit("The electric field contains nonfinite values")

    fig = plot_ez(space.group_nodes, ez_groups)
    if output is None:
        plt.show()
    else:
        fig.savefig(output)
        print(f"Saved plot to {output}", flush=True)

    if max_error >= 0.03:
        raise SystemExit("The electric field exceeded the cavity error tolerance")


if __name__ == "__main__":
    parser = ArgumentParser(description="Plot the electric field in a 2D PEC cavity.")
    parser.add_argument("--output", help="Save the plot to a file instead of opening a window")
    main(parser.parse_args().output)
