# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Enable coverage flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

# Check for the gcovr version required by the coverage stylesheet.
find_program(GCOVR_EXECUTABLE gcovr)
set(OISEAU_GCOVR_VERSION_OK FALSE)
if(GCOVR_EXECUTABLE)
    execute_process(
        COMMAND "${GCOVR_EXECUTABLE}" --version
        OUTPUT_VARIABLE GCOVR_VERSION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE GCOVR_VERSION_RESULT
    )
    string(REGEX MATCH "^gcovr [0-9]+[.][0-9]+" GCOVR_VERSION_LINE "${GCOVR_VERSION_OUTPUT}")
    if(GCOVR_VERSION_RESULT EQUAL 0 AND "${GCOVR_VERSION_LINE}" STREQUAL "gcovr 7.2")
        set(OISEAU_GCOVR_VERSION_OK TRUE)
    endif()
endif()

if(NOT OISEAU_GCOVR_VERSION_OK)
    if(NOT OISEAU_COVERAGE_AUTO_INSTALL_GCOVR)
        message(FATAL_ERROR "Coverage reporting requires gcovr 7.2; found: ${GCOVR_VERSION_LINE}. "
                            "Install gcovr 7.2 or enable OISEAU_COVERAGE_AUTO_INSTALL_GCOVR."
        )
    endif()

    find_package(
        Python3
        COMPONENTS Interpreter
        REQUIRED
    )
    set(OISEAU_GCOVR_VENV "${CMAKE_BINARY_DIR}/_deps/gcovr-7.2")
    set(OISEAU_GCOVR_PYTHON "${OISEAU_GCOVR_VENV}/bin/python")
    set(OISEAU_GCOVR_EXECUTABLE "${OISEAU_GCOVR_VENV}/bin/gcovr")
    set(OISEAU_GCOVR_VENV_OK FALSE)

    if(EXISTS "${OISEAU_GCOVR_EXECUTABLE}")
        execute_process(
            COMMAND "${OISEAU_GCOVR_EXECUTABLE}" --version
            OUTPUT_VARIABLE OISEAU_GCOVR_VENV_VERSION_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE OISEAU_GCOVR_VENV_VERSION_RESULT
        )
        string(REGEX MATCH "^gcovr [0-9]+[.][0-9]+" OISEAU_GCOVR_VENV_VERSION_LINE
                     "${OISEAU_GCOVR_VENV_VERSION_OUTPUT}"
        )
        if(OISEAU_GCOVR_VENV_VERSION_RESULT EQUAL 0 AND "${OISEAU_GCOVR_VENV_VERSION_LINE}"
                                                        STREQUAL "gcovr 7.2"
        )
            set(OISEAU_GCOVR_VENV_OK TRUE)
        endif()
    endif()

    if(NOT OISEAU_GCOVR_VENV_OK)
        message(STATUS "Installing gcovr 7.2 into ${OISEAU_GCOVR_VENV}")
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -m venv "${OISEAU_GCOVR_VENV}"
            RESULT_VARIABLE OISEAU_GCOVR_VENV_RESULT
            ERROR_VARIABLE OISEAU_GCOVR_VENV_ERROR
        )
        if(NOT OISEAU_GCOVR_VENV_RESULT EQUAL 0)
            message(
                FATAL_ERROR
                    "Could not create the gcovr virtual environment: ${OISEAU_GCOVR_VENV_ERROR}"
            )
        endif()

        execute_process(
            COMMAND "${OISEAU_GCOVR_PYTHON}" -m pip install --disable-pip-version-check --upgrade
                    "gcovr==7.2"
            RESULT_VARIABLE OISEAU_GCOVR_INSTALL_RESULT
            ERROR_VARIABLE OISEAU_GCOVR_INSTALL_ERROR
        )
        if(NOT OISEAU_GCOVR_INSTALL_RESULT EQUAL 0)
            message(
                FATAL_ERROR
                    "Could not install gcovr 7.2: ${OISEAU_GCOVR_INSTALL_ERROR}. "
                    "Check Python's venv/pip support and network access, or install gcovr 7.2 manually."
            )
        endif()
    endif()

    set(GCOVR_EXECUTABLE "${OISEAU_GCOVR_EXECUTABLE}")
    execute_process(
        COMMAND "${GCOVR_EXECUTABLE}" --version
        OUTPUT_VARIABLE GCOVR_VERSION_OUTPUT
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE GCOVR_VERSION_RESULT
    )
    string(REGEX MATCH "^gcovr [0-9]+[.][0-9]+" GCOVR_VERSION_LINE "${GCOVR_VERSION_OUTPUT}")
    if(NOT GCOVR_VERSION_RESULT EQUAL 0 OR NOT "${GCOVR_VERSION_LINE}" STREQUAL "gcovr 7.2")
        message(FATAL_ERROR "Expected gcovr 7.2 after installation; found: ${GCOVR_VERSION_LINE}")
    endif()
endif()

message(STATUS "gcovr found: ${GCOVR_EXECUTABLE} (${GCOVR_VERSION_LINE})")

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
    COMMAND "${GCOVR_EXECUTABLE}" --delete --filter ${CMAKE_SOURCE_DIR}/src/ --html-details -o
            ${COVERAGE_DIR}/index.html
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/website/coverage/custom.css
            ${COVERAGE_DIR}/index.css
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
