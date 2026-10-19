# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)
project(oiseau_pybind11)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(OISEAU_BUILD_DEMO OFF)
set(OISEAU_BUILD_TEST OFF)
set(OISEAU_BUILD_DOCS OFF)
set(OISEAU_BUILD_BENCHMARK OFF)
set(OISEAU_BUILD_SHARED OFF)

add_subdirectory(.. oiseau)
pybind11_add_module(
    cpp MODULE oiseau/wrappers/oiseau.cpp oiseau/wrappers/mesh.cpp oiseau/wrappers/io.cpp
)
target_link_libraries(cpp PRIVATE oiseau)
target_include_directories(cpp PRIVATE ../src)
