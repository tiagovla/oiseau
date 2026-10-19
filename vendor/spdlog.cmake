# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG ${SPDLOG_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
)
FetchContent_MakeAvailable(spdlog)

set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
