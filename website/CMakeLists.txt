# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(WEBSITE_DIR ${CMAKE_BINARY_DIR}/web)
set(DOCS_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs/html)
set(COVERAGE_OUTPUT_DIR ${CMAKE_BINARY_DIR}/coverage)

add_custom_target(
    website
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${WEBSITE_DIR}
    COMMAND npm run build -- --outDir ${WEBSITE_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOCS_OUTPUT_DIR} ${WEBSITE_DIR}/docs
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building the Astro website and copying documentation"
    VERBATIM
)

add_dependencies(website docs)
if(TARGET coverage)
    add_dependencies(website coverage)
    add_custom_command(
        TARGET website
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${COVERAGE_OUTPUT_DIR} ${WEBSITE_DIR}/coverage
        COMMENT "Copying coverage files"
        VERBATIM
    )
endif()

add_custom_command(
    TARGET website
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "To serve the website,"
    COMMAND ${CMAKE_COMMAND} -E echo "run python -m http.server 8080 --directory ${WEBSITE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "and open http://localhost:8080 in your browser"
    COMMENT "Serving the website"
)
