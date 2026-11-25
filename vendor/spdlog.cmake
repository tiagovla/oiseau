# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

set(OISEAU_LOG_LEVEL
    "INFO"
    CACHE STRING "Default runtime logging level"
)
set(OISEAU_LOG_LEVELS
    TRACE
    DEBUG
    INFO
    WARN
    ERROR
    CRITICAL
    OFF
)
set_property(CACHE OISEAU_LOG_LEVEL PROPERTY STRINGS ${OISEAU_LOG_LEVELS})
string(TOUPPER "${OISEAU_LOG_LEVEL}" OISEAU_LOG_LEVEL)
if(NOT OISEAU_LOG_LEVEL IN_LIST OISEAU_LOG_LEVELS)
    message(
        FATAL_ERROR
            "Invalid OISEAU_LOG_LEVEL '${OISEAU_LOG_LEVEL}'. Choose one of ${OISEAU_LOG_LEVELS}."
    )
endif()
set(OISEAU_LOG_LEVEL
    "${OISEAU_LOG_LEVEL}"
    PARENT_SCOPE
)

set(SPDLOG_FMT_EXTERNAL ON)

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG ${SPDLOG_VERSION}
    GIT_PROGRESS 1
)
FetchContent_MakeAvailable(spdlog)

set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
