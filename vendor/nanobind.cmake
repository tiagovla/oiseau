# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

if(OISEAU_BUILD_PYTHON_EXTENSION)
    find_package(
        Python 3.13 REQUIRED COMPONENTS Interpreter Development.Module Development.SABIModule
    )
else()
    find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development)
    message(STATUS "Python environment for C++ plotting: ${Python_EXECUTABLE}")

    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import matplotlib, PyQt6"
        RESULT_VARIABLE PLOTTING_IMPORT_RESULT
        ERROR_QUIET
    )
    if(NOT PLOTTING_IMPORT_RESULT EQUAL 0)
        message(FATAL_ERROR "Matplotlib and PyQt6 are required by ${Python_EXECUTABLE}")
    endif()
endif()

include(FetchContent)
FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v3.1.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(nanobind)

if(NOT OISEAU_BUILD_PYTHON_EXTENSION)
    set(OISEAU_NANOBIND_TARGET nanobind-static)
    if(OISEAU_BUILD_SHARED)
        # The library and its callers must use the same nanobind runtime.
        set(OISEAU_NANOBIND_TARGET nanobind)
        if(APPLE)
            set_property(
                TARGET oiseau
                APPEND
                PROPERTY INSTALL_RPATH "@loader_path"
            )
        elseif(UNIX)
            set_property(
                TARGET oiseau
                APPEND
                PROPERTY INSTALL_RPATH "$ORIGIN"
            )
        endif()
    endif()
    nanobind_build_library(${OISEAU_NANOBIND_TARGET})
    target_link_libraries(${OISEAU_NANOBIND_TARGET} PRIVATE Python::Python)
    target_link_libraries(oiseau_deps INTERFACE ${OISEAU_NANOBIND_TARGET} Python::Python)
    install(
        TARGETS ${OISEAU_NANOBIND_TARGET}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
    )
endif()
