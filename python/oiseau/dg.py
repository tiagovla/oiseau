# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Discontinuous Galerkin spaces, equations, and time stepping."""

from .cpp.dg import (
    DGSpace,
    DGState,
    Equation,
    EquationSystem,
    Field,
    LinearDG2D,
    LinearExpression,
    dt,
    dx,
    dy,
    run_rk4,
)

__all__ = [
    "DGSpace",
    "DGState",
    "Equation",
    "EquationSystem",
    "Field",
    "LinearDG2D",
    "LinearExpression",
    "dt",
    "dx",
    "dy",
    "run_rk4",
]
