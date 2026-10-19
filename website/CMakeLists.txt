# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(WEBSITE_DIR ${CMAKE_BINARY_DIR}/web)
set(SOURCE_WEBSITE_DIR ${CMAKE_SOURCE_DIR}/website/oiseau)
set(DOCS_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs/html)
set(COVERAGE_OUTPUT_DIR ${CMAKE_BINARY_DIR}/coverage)

add_custom_target(website COMMENT "Generating website content")

add_custom_command(
    TARGET website
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${WEBSITE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${WEBSITE_DIR}
    COMMENT "Cleaning and creating website directory"
)

add_dependencies(website docs coverage)

add_custom_command(
    TARGET website
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${SOURCE_WEBSITE_DIR} ${WEBSITE_DIR}
    COMMENT "Copying website source files"
)
#
add_custom_command(
    TARGET website
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${COVERAGE_OUTPUT_DIR} ${WEBSITE_DIR}/coverage
    COMMENT "Copying coverage files"
)

add_custom_command(
    TARGET website
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOCS_OUTPUT_DIR} ${WEBSITE_DIR}/docs
    COMMENT "Copying documentation files"
)

add_custom_command(
    TARGET website
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "To serve the website,"
    COMMAND ${CMAKE_COMMAND} -E echo "run python -m http.server 8080 --directory ${WEBSITE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "and open http://localhost:8080 in your browser"
    COMMENT "Serving the website"
)
