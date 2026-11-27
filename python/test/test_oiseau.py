# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
import pytest

from oiseau.dg import EquationSystem, Field, LinearDG2D, Space, dt, dx, run_rk4
from oiseau.mesh import CellKind, Mesh


@pytest.mark.parametrize("boundary_state", [None, lambda interior, normal: interior])
def test_linear_advection_preserves_a_constant_state(boundary_state):
    mesh = Mesh([(0, 0), (1, 0), (0, 1)], [[0, 1, 2]], [CellKind.Triangle])
    space = Space(mesh, order=1)
    field = Field("u")
    equations = EquationSystem([field], [dt(field) == -dx(field)])
    solver = LinearDG2D(mesh, space, equations, wave_speed=1.0, boundary_state=boundary_state)
    initial = solver.initial_state(field, lambda x, y: 2.0)

    (values,) = run_rk4(solver, initial, initial_time=0.0, time_step=0.01, steps=1).groups

    np.testing.assert_allclose(values, np.full((1, 3, 1), 2.0), rtol=0, atol=1e-12)
