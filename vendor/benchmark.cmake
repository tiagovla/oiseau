# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(BENCHMARK_ENABLE_TESTING OFF)

FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG ${BENCHMARK_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
)

FetchContent_MakeAvailable(benchmark)
