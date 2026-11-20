# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

include(CheckFunctionExists)

function(oiseau_find_system_blas)
    if(CMAKE_IMPORT_LIBRARY_SUFFIX)
        set(CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_IMPORT_LIBRARY_SUFFIX}")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    set(CMAKE_REQUIRED_QUIET TRUE)

    find_library(OISEAU_SHARED_OPENBLAS NAMES openblas NO_CACHE)
    if(OISEAU_SHARED_OPENBLAS)
        set(CMAKE_REQUIRED_LIBRARIES "${OISEAU_SHARED_OPENBLAS}")
        unset(OISEAU_OPENBLAS_HAS_LAPACK CACHE)
        unset(OISEAU_OPENBLAS_HAS_CBLAS CACHE)
        check_function_exists(dgeev_ OISEAU_OPENBLAS_HAS_LAPACK)
        check_function_exists(cblas_dgemm OISEAU_OPENBLAS_HAS_CBLAS)
        if(OISEAU_OPENBLAS_HAS_LAPACK AND OISEAU_OPENBLAS_HAS_CBLAS)
            target_link_libraries(xtensor_stack INTERFACE "${OISEAU_SHARED_OPENBLAS}")
            set(OISEAU_BLAS_PROVIDER
                "shared OpenBLAS"
                PARENT_SCOPE
            )
            return()
        endif()
    endif()

    set(CMAKE_REQUIRED_LIBRARIES "")
    set(BLA_STATIC OFF)
    set(BLA_SIZEOF_INTEGER 4)
    find_package(LAPACK QUIET MODULE)
    if(NOT LAPACK_FOUND)
        return()
    endif()

    set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
    unset(OISEAU_SYSTEM_HAS_CBLAS CACHE)
    check_function_exists(cblas_dgemm OISEAU_SYSTEM_HAS_CBLAS)
    if(NOT OISEAU_SYSTEM_HAS_CBLAS)
        find_library(OISEAU_SYSTEM_CBLAS NAMES cblas NO_CACHE)
        if(NOT OISEAU_SYSTEM_CBLAS)
            return()
        endif()
        list(PREPEND CMAKE_REQUIRED_LIBRARIES "${OISEAU_SYSTEM_CBLAS}")
        unset(OISEAU_SEPARATE_CBLAS_WORKS CACHE)
        check_function_exists(cblas_dgemm OISEAU_SEPARATE_CBLAS_WORKS)
        if(NOT OISEAU_SEPARATE_CBLAS_WORKS)
            return()
        endif()
    endif()

    target_link_libraries(xtensor_stack INTERFACE ${CMAKE_REQUIRED_LIBRARIES})
    set(OISEAU_BLAS_PROVIDER
        "system BLAS/LAPACK"
        PARENT_SCOPE
    )
endfunction()

oiseau_find_system_blas()
if(NOT OISEAU_BLAS_PROVIDER)
    include(openblas.cmake)
    target_link_libraries(xtensor_stack INTERFACE openblas_static)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        # Keep BLAS symbols private when loading Oiseau alongside NumPy.
        target_link_options(xtensor_stack INTERFACE "LINKER:--exclude-libs,libopenblas.a")
    endif()
    set(OISEAU_BLAS_PROVIDER "bundled static OpenBLAS")
endif()
message(STATUS "Oiseau BLAS/LAPACK: ${OISEAU_BLAS_PROVIDER}")
