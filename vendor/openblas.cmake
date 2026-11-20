# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

function(oiseau_add_openblas)
    set(BUILD_STATIC_LIBS ON)
    set(BUILD_SHARED_LIBS OFF)
    set(BUILD_TESTING OFF)
    set(BUILD_WITHOUT_LAPACK OFF)
    set(BUILD_WITHOUT_LAPACKE ON)
    set(C_LAPACK ON)
    set(NOFORTRAN ON)
    set(DYNAMIC_ARCH ON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    FetchContent_Declare(
        openblas
        GIT_REPOSITORY https://github.com/OpenMathLib/OpenBLAS.git
        GIT_TAG ${OPENBLAS_VERSION}
        GIT_SHALLOW 1
    )
    FetchContent_MakeAvailable(openblas)
endfunction()

oiseau_add_openblas()
