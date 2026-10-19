// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/utils/logging.hpp"

#include <spdlog/pattern_formatter.h>

#include <memory>

#include "spdlog/cfg/argv.h"
#include "spdlog/common.h"
#include "spdlog/logger.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

namespace oiseau::logging {

static std::shared_ptr<spdlog::logger> default_logger;

void set_verbosity(Verbosity level) {
  spdlog::set_level(static_cast<spdlog::level::level_enum>(level));
}

void init(int argc, char** argv) {
  spdlog::cfg::load_argv_levels(argc, argv);

  spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%# %!] %v");
  if (!default_logger) {
    default_logger = spdlog::stdout_color_mt("oiseau");
    spdlog::set_default_logger(default_logger);
  }
  spdlog::set_level(spdlog::level::info);
}

}  // namespace oiseau::logging
