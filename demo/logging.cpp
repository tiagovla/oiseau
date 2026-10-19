// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/utils/logging.hpp"

#include <spdlog/cfg/argv.h>

#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

int main(int argc, char** argv) {
  oiseau::logging::init(argc, argv);
  oiseau::logging::set_verbosity(oiseau::logging::Verbosity::DEBUG);

  xt::xarray<double> arr = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  LOG_DEBUG("Array shape: {} x {}", arr.shape()[0], arr.shape()[1]);
  LOG_DEBUG("{}", arr);
}
