# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

macro(oiseau_add_executable name)
    add_executable(oiseau_${name} ${name}.cpp)
    target_link_libraries(oiseau_${name} PRIVATE oiseau oiseau_deps)
    set_target_properties(
        oiseau_${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/demo
    )
endmacro()

oiseau_add_executable(nodaldg)
oiseau_add_executable(dg_connectivity)
oiseau_add_executable(dg_space)
oiseau_add_executable(mesh)
oiseau_add_executable(io)
oiseau_add_executable(plotting)
oiseau_add_executable(logging)
