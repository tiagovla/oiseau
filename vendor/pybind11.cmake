# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(PYBIND11_FINDPYTHON "NEW")

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG ${PYBIND11_VERSION}
    GIT_SHALLOW 0
    GIT_PROGRESS 1
)

FetchContent_MakeAvailable(pybind11)

message(STATUS "Python environment found by pybind11: ${Python_EXECUTABLE}")

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import matplotlib"
    RESULT_VARIABLE MATPLOTLIB_RESULT
    ERROR_QUIET
)
if(NOT MATPLOTLIB_RESULT EQUAL 0)
    message(FATAL_ERROR "Matplotlib is required by ${Python_EXECUTABLE}")
endif()

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import PyQt6"
    RESULT_VARIABLE PYQT6_RESULT
    ERROR_QUIET
)
if(NOT PYQT6_RESULT EQUAL 0)
    message(FATAL_ERROR "PyQt6 is required by ${Python_EXECUTABLE}")
endif()
