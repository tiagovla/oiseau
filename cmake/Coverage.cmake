# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Enable coverage flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

# Check for gcovr
find_program(GCOVR_EXECUTABLE gcovr)
if(NOT GCOVR_EXECUTABLE)
    message(FATAL_ERROR "gcovr not found! Please install it to enable coverage reporting.")
else()
    message(STATUS "gcovr found: ${GCOVR_EXECUTABLE}")
endif()

# Coverage output directory
set(COVERAGE_DIR ${CMAKE_BINARY_DIR}/coverage)

# Define the custom target
add_custom_target(coverage COMMENT "Building and generating coverage report")

# PRE_BUILD: clean previous coverage data
add_custom_command(
    TARGET coverage
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning previous coverage data..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${COVERAGE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${COVERAGE_DIR}
)

# POST_BUILD: run tests and generate report
add_custom_command(
    TARGET coverage
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Running tests and generating coverage report..."
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
    COMMAND ${GCOVR_EXECUTABLE} --delete --filter ${CMAKE_SOURCE_DIR}/src/ --html-details -o
            ${COVERAGE_DIR}/index.html
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/website/coverage/custom.css
            ${COVERAGE_DIR}/index.css
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
