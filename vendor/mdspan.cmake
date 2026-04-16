# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

FetchContent_Declare(
    mdspan
    GIT_REPOSITORY https://github.com/kokkos/mdspan.git
    GIT_TAG ${MDSPAN_VERSION}
    GIT_SHALLOW 1
    GIT_PROGRESS 1
)

FetchContent_MakeAvailable(mdspan)

add_library(std::mdspan ALIAS mdspan)
# target_compile_definitions(mdspan INTERFACE MDSPAN_IMPL_STANDARD_NAMESPACE=std)
# target_compile_definitions(mdspan INTERFACE MDSPAN_IMPL_PROPOSED_NAMESPACE=experimental)
