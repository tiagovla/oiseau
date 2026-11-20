# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.26)
project(oiseau_python)

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

include(../vendor/nanobind.cmake)

nanobind_add_module(
    cpp
    NB_STATIC
    STABLE_ABI
    oiseau/wrappers/oiseau.cpp
    oiseau/wrappers/mesh.cpp
    oiseau/wrappers/io.cpp
    oiseau/wrappers/dg.cpp
)
if(NOT TARGET nanobind-static-abi3)
    message(FATAL_ERROR "The stable ABI build requires standard CPython 3.13 or later")
endif()
if(NOT WIN32)
    # Use the stable ABI suffix that Python 3.13 and 3.14 can load.
    set_target_properties(cpp PROPERTIES SUFFIX ".abi3${CMAKE_SHARED_MODULE_SUFFIX}")
endif()
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
