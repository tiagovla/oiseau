# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_test(oiseau_test_dg_nodal_utils test_utils.cpp)

add_test(oiseau_test_dg_nodal_ref_element test_ref_element.cpp)
add_test(oiseau_test_dg_nodal_ref_line test_ref_line.cpp)
add_test(oiseau_test_dg_nodal_ref_triangle test_ref_triangle.cpp)
add_test(oiseau_test_dg_nodal_ref_quadrilateral test_ref_quadrilateral.cpp)
add_test(oiseau_test_dg_nodal_ref_tetrahedron test_ref_tetrahedron.cpp)
add_test(oiseau_test_dg_nodal_ref_hexahedron test_ref_hexahedron.cpp)
