# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

macro(add_test test_name)
    add_executable(${test_name} ${ARGN})
    target_link_libraries(${test_name} PRIVATE gtest gtest_main oiseau oiseau_deps)
    set_target_properties(
        ${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test
    )
    target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/test)
    gtest_discover_tests(${test_name})
endmacro()
