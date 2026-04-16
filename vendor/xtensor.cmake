# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

find_package(LAPACK REQUIRED)
find_package(OpenMP QUIET)
find_package(OpenBLAS QUIET)

FetchContent_Declare(
    xtl
    GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
    GIT_TAG ${XTL_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
    FIND_PACKAGE_ARGS 0.8.0
)

FetchContent_Declare(
    xtensor
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
    GIT_TAG ${XTENSOR_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
    FIND_PACKAGE_ARGS 0.26.0
)

FetchContent_Declare(
    xtensor-blas
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-blas
    GIT_TAG ${XTENSOR-BLAS_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
    FIND_PACKAGE_ARGS 0.22.0
)

FetchContent_MakeAvailable(xtl xtensor xtensor-blas)

add_library(xtensor_stack INTERFACE)
target_link_libraries(xtensor_stack INTERFACE xtl xtensor xtensor-blas LAPACK::LAPACK)

if(OpenBLAS_FOUND AND OpenMP_FOUND)
    target_compile_definitions(xtensor-blas INTERFACE WITH_OPENBLAS)
    target_link_libraries(xtensor_stack INTERFACE OpenBLAS::OpenBLAS)
    message(STATUS "xtensor_stack: OpenBLAS found, linking against OpenBLAS")
endif()

get_target_property(XTL_INCLUDES xtl INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(XTENSOR_INCLUDES xtensor INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(XTENSOR_BLAS_INCLUDES xtensor-blas INTERFACE_INCLUDE_DIRECTORIES)

target_include_directories(
    xtensor_stack SYSTEM INTERFACE ${XTL_INCLUDES} ${XTENSOR_INCLUDES} ${XTENSOR_BLAS_INCLUDES}
)
