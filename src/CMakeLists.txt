# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

file(
    GLOB_RECURSE
    OISEAU_CORE_SOURCES
    CONFIGURE_DEPENDS
    "oiseau/dg/*.cpp"
    "oiseau/io/*.cpp"
    "oiseau/mesh/*.cpp"
    "oiseau/utils/*.cpp"
)
if(NOT OISEAU_BUILD_PYTHON_EXTENSION)
    file(GLOB OISEAU_PLOTTING_SOURCES CONFIGURE_DEPENDS "oiseau/plotting/*.cpp")
    list(APPEND OISEAU_CORE_SOURCES ${OISEAU_PLOTTING_SOURCES})
endif()

set(OISEAU_PUBLIC_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}
    PARENT_SCOPE
)
set(OISEAU_CORE_SOURCES
    ${OISEAU_CORE_SOURCES}
    PARENT_SCOPE
)
