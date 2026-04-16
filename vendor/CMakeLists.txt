# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

include(FetchContent)

# Library version pins
set(BENCHMARK_VERSION v1.9.3)
set(BLAZE_VERSION v3.8.2)
set(BOOST_VERSION boost-1.88.0)
set(EIGEN_VERSION 171bd08ca987987c3c50f0fa5dd8914bdd42dd3b)
set(FMT_VERSION 11.2.0)
set(GMSH_VERSION gmsh_4_13_1)
set(GOOGLETEST_VERSION v1.17.0)
set(MDSPAN_VERSION stable)
set(PYBIND11_VERSION 67424358f457a4957785f6b0912b73fdec962ec9)
set(SPDLOG_VERSION v1.15.3)
set(XTENSOR-BLAS_VERSION 0.22.0)
set(XTENSOR_VERSION 0.26.0)
set(XTL_VERSION 0.8.0)

# FetchContent options
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED OFF)

# Optional dependencies
if(OISEAU_BUILD_TEST)
    include(googletest.cmake)
endif()

if(OISEAU_BUILD_BENCHMARK)
    include(benchmark.cmake)
endif()

if(OISEAU_BUILD_DOCS)
    include(doxygen.cmake)
endif()

# Create interface target for dependencies
add_library(oiseau_deps INTERFACE)

# Core dependencies
include(xtensor.cmake)
include(fmt.cmake)
include(spdlog.cmake)
include(pybind11.cmake)
include(mdspan.cmake)

target_link_libraries(
    oiseau_deps INTERFACE xtensor_stack fmt::fmt spdlog::spdlog pybind11::embed std::mdspan
)
