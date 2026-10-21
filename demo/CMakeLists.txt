# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

file(GLOB OISEAU_DEMO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(demo_source IN LISTS OISEAU_DEMO_SOURCES)
    get_filename_component(demo_name "${demo_source}" NAME_WE)
    # if(demo_name MATCHES "^dev") continue() endif()
    add_executable(oiseau_${demo_name} "${demo_source}")
    target_link_libraries(oiseau_${demo_name} PRIVATE oiseau oiseau_deps)
    set_target_properties(
        oiseau_${demo_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/demo
    )
endforeach()
