# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(MAKEINFO OFF)
set(ENABLE_PRIVATE_API ON)

string(REPLACE "." "_" gmsh_tag_version "${GMSH_VERSION}")

FetchContent_Declare(
    gmsh
    GIT_REPOSITORY https://gitlab.onelab.info/gmsh/gmsh
    GIT_TAG gmsh_${gmsh_tag_version}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
    FIND_PACKAGE_ARGS ${GMSH_VERSION}
)

FetchContent_MakeAvailable(gmsh)
