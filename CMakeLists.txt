# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.14)

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
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ------------------------------------------------------------------------------
# Build Options
# ------------------------------------------------------------------------------

option(OISEAU_BUILD_DEMO "Build demo applications" ON)
option(OISEAU_BUILD_TEST "Build unit tests" ON)
option(OISEAU_BUILD_DOCS "Build documentation" OFF)
option(OISEAU_BUILD_BENCHMARK "Build benchmarks" OFF)
option(OISEAU_BUILD_COVERAGE "Build with coverage" OFF)
option(OISEAU_BUILD_SHARED "Build shared library" OFF)

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
    message(STATUS "oiseau: Building SHARED library")
    add_library(oiseau SHARED ${OISEAU_CORE_SOURCES})
else()
    message(STATUS "oiseau: Building STATIC library")
    add_library(oiseau STATIC ${OISEAU_CORE_SOURCES})
    target_compile_options(oiseau PRIVATE -fPIC)
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
    message(STATUS "oiseau: Coverage build enabled")
    include(Coverage)
endif()

if(OISEAU_BUILD_DEMO)
    message(STATUS "oiseau: Demo build enabled")
    add_subdirectory(demo)
endif()

if(OISEAU_BUILD_TEST)
    message(STATUS "oiseau: Test build enabled")
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()

if(OISEAU_BUILD_BENCHMARK)
    message(STATUS "oiseau: Benchmark build enabled")
    add_subdirectory(benchmark)
endif()

if(OISEAU_BUILD_DOCS)
    message(STATUS "oiseau: Documentation build enabled")
    add_subdirectory(doc)
endif()

# ------------------------------------------------------------------------------
# Website
# ------------------------------------------------------------------------------
if(OISEAU_BUILD_DOCS AND OISEAU_BUILD_TEST)
    add_subdirectory(website)
endif()
