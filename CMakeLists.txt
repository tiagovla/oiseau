# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.16)

project(
    oiseau
    VERSION 0.1.0
    LANGUAGES C CXX
)

# ------------------------------------------------------------------------------
# Build Configuration
# ------------------------------------------------------------------------------

# Enable compile commands for tools like clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require modern C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Debug build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        "Debug"
        CACHE STRING "Choose the type of build: Debug, Release, RelWithDebInfo, MinSizeRel." FORCE
    )
endif()

# ------------------------------------------------------------------------------
# Build Options
# ------------------------------------------------------------------------------

option(OISEAU_BUILD_DEMO "Build demo applications" ON)
option(OISEAU_BUILD_TEST "Build unit tests" ON)
option(OISEAU_BUILD_DOCS "Build documentation" OFF)
option(OISEAU_BUILD_BENCHMARK "Build benchmarks" OFF)
option(OISEAU_BUILD_COVERAGE "Build with coverage" OFF)
option(OISEAU_COVERAGE_AUTO_INSTALL_GCOVR
       "Install the required gcovr version into the build directory when needed" ON
)
option(OISEAU_BUILD_SHARED "Build shared library" OFF)
option(OISEAU_PRECOMPILE_HEADERS "Precompile common STL and xtensor headers" OFF)

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

add_subdirectory(vendor EXCLUDE_FROM_ALL)

# ------------------------------------------------------------------------------
# Library Target
# ------------------------------------------------------------------------------

add_subdirectory(src)

# conditional to build static or shared
if(OISEAU_BUILD_SHARED)
    add_library(oiseau SHARED ${OISEAU_CORE_SOURCES})
else()
    add_library(oiseau STATIC ${OISEAU_CORE_SOURCES})
    target_compile_options(oiseau PRIVATE -fPIC)
endif()

if(OISEAU_PRECOMPILE_HEADERS)
    target_precompile_headers(
        oiseau
        PRIVATE
        <algorithm>
        <array>
        <cmath>
        <cstddef>
        <cstdint>
        <limits>
        <memory>
        <span>
        <stdexcept>
        <string>
        <string_view>
        <utility>
        <vector>
        <xtensor-blas/xlinalg.hpp>
        <xtensor/containers/xarray.hpp>
        <xtensor/containers/xtensor.hpp>
        <xtensor/core/xtensor_forward.hpp>
        <xtensor/generators/xbuilder.hpp>
        <xtensor/misc/xmanipulation.hpp>
        <xtensor/views/xslice.hpp>
        <xtensor/views/xview.hpp>
    )
endif()

target_link_libraries(oiseau PRIVATE oiseau_deps)
target_include_directories(
    oiseau PUBLIC $<BUILD_INTERFACE:${OISEAU_PUBLIC_INCLUDE_DIR}> $<INSTALL_INTERFACE:include>
)

# ------------------------------------------------------------------------------
# Install the library
# ------------------------------------------------------------------------------

install(
    TARGETS oiseau
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES
    DESTINATION include
)
install(
    DIRECTORY ${OISEAU_PUBLIC_INCLUDE_DIR}/
    DESTINATION include
    FILES_MATCHING
    PATTERN "*.hpp"
)

# CMAKE Modules

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ------------------------------------------------------------------------------
# Optional Subdirectories
# ------------------------------------------------------------------------------
if(OISEAU_BUILD_COVERAGE)
    include(Coverage)
endif()

if(OISEAU_BUILD_DEMO)
    add_subdirectory(demo)
endif()

if(OISEAU_BUILD_TEST)
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()

if(OISEAU_BUILD_BENCHMARK)
    add_subdirectory(benchmark)
endif()

if(OISEAU_BUILD_DOCS)
    add_subdirectory(doc)
endif()

# ------------------------------------------------------------------------------
# Website
# ------------------------------------------------------------------------------
if(OISEAU_BUILD_DOCS AND OISEAU_BUILD_TEST)
    add_subdirectory(website)
endif()

set(OISEAU_ENABLED_FEATURES "")
if(OISEAU_BUILD_DEMO)
    list(APPEND OISEAU_ENABLED_FEATURES demos)
endif()
if(OISEAU_BUILD_TEST)
    list(APPEND OISEAU_ENABLED_FEATURES tests)
endif()
if(OISEAU_BUILD_BENCHMARK)
    list(APPEND OISEAU_ENABLED_FEATURES benchmarks)
endif()
if(OISEAU_BUILD_DOCS)
    list(APPEND OISEAU_ENABLED_FEATURES docs)
endif()
if(OISEAU_BUILD_COVERAGE)
    list(APPEND OISEAU_ENABLED_FEATURES coverage)
endif()
if(OISEAU_ENABLED_FEATURES)
    list(JOIN OISEAU_ENABLED_FEATURES ", " OISEAU_FEATURES_TEXT)
else()
    set(OISEAU_FEATURES_TEXT "core only")
endif()

if(OISEAU_BUILD_SHARED)
    set(OISEAU_LIBRARY_KIND "shared")
else()
    set(OISEAU_LIBRARY_KIND "static")
endif()

string(CONCAT OISEAU_CONFIG_SUMMARY "  version   ${PROJECT_VERSION}\n"
              "  build     ${CMAKE_BUILD_TYPE} / ${OISEAU_LIBRARY_KIND}\n"
              "  enabled   ${OISEAU_FEATURES_TEXT}"
)
if(OISEAU_BUILD_COVERAGE)
    string(APPEND OISEAU_CONFIG_SUMMARY "\n  coverage  ${GCOVR_VERSION_LINE}")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/OiseauBanner.cmake")
message(NOTICE "${OISEAU_CONFIG_SUMMARY}")
