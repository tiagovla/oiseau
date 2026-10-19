# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

macro(add_benchmark benchmark_name)
    add_executable(${benchmark_name} ${ARGN})
    target_link_libraries(${benchmark_name} PRIVATE oiseau oiseau_deps)
    target_link_libraries(${benchmark_name} PRIVATE benchmark::benchmark benchmark::benchmark_main)
    set_target_properties(
        ${benchmark_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks
    )
endmacro()

add_benchmark(oiseau_benchmark_oiseau benchmark_oiseau.cpp)
add_benchmark(oiseau_benchmark_xtensor benchmark_xtensor.cpp)
add_benchmark(oiseau_benchmark_dot_layout benchmark_dot_layout.cpp)
