# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG ${EIGEN_VERSION}
    GIT_SHALLOW 0
    GIT_PROGRESS 1
)

FetchContent_MakeAvailable(Eigen3)
