# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

FetchContent_Declare(
    Boost
    GIT_REPOSITORY https://github.com/boostorg/boost
    GIT_TAG ${BOOST_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
)

FetchContent_MakeAvailable(Boost)
