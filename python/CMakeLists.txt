# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)
project(oiseau_pybind11)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(OISEAU_BUILD_DEMO OFF)
set(OISEAU_BUILD_TEST OFF)
set(OISEAU_BUILD_DOCS OFF)
set(OISEAU_BUILD_BENCHMARK OFF)
set(OISEAU_BUILD_SHARED OFF)
set(OISEAU_BUILD_PYTHON_EXTENSION ON)

add_subdirectory(.. oiseau)
pybind11_add_module(
    cpp MODULE oiseau/wrappers/oiseau.cpp oiseau/wrappers/mesh.cpp oiseau/wrappers/io.cpp
    oiseau/wrappers/dg.cpp
)
target_link_libraries(cpp PRIVATE oiseau xtensor_stack)
target_include_directories(cpp PRIVATE ../src)

if(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(WHEEL_LICENSE_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/licenses")
    configure_file(../LICENSE "${WHEEL_LICENSE_DIR}/Oiseau.txt" COPYONLY)
    if(TARGET openblas_static)
        get_target_property(OPENBLAS_SOURCE_DIR openblas_static SOURCE_DIR)
        configure_file(
            "${OPENBLAS_SOURCE_DIR}/LICENSE" "${WHEEL_LICENSE_DIR}/OpenBLAS.txt" COPYONLY
        )
        configure_file(
            "${OPENBLAS_SOURCE_DIR}/lapack-netlib/LICENSE" "${WHEEL_LICENSE_DIR}/LAPACK.txt"
            COPYONLY
        )
    endif()
endif()
