# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_test(oiseau_test_dg_geometry test_dg_geometry.cpp)
add_test(oiseau_test_dg_field test_field.cpp)
add_test(oiseau_test_dg_equation test_equation.cpp)

add_subdirectory(nodal)
