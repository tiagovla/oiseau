// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#if defined(__cpp_lib_mdspan)
#include <mdspan>
#else
#include <mdspan/mdspan.hpp>

#include "experimental/__p0009_bits/layout_stride.hpp"
namespace std {
using Kokkos::dextents;
using Kokkos::extents;
using Kokkos::layout_right;
using Kokkos::mdspan;
}  // namespace std
#endif

#include <experimental/__p1684_bits/mdarray.hpp>
namespace stdex {
using Kokkos::Experimental::mdarray;

}  // namespace stdex
