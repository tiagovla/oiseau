// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <spdlog/cfg/argv.h>
#include <spdlog/common.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <format>
#include <sstream>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/io/xio.hpp>

#undef SPDLOG_TRACE
#undef SPDLOG_DEBUG
#undef SPDLOG_INFO
#undef SPDLOG_WARN
#undef SPDLOG_ERROR
#undef SPDLOG_CRITICAL

namespace oiseau::logging {

enum class Verbosity {
  OFF = SPDLOG_LEVEL_OFF,
  FATAL = SPDLOG_LEVEL_CRITICAL,
  ERROR = SPDLOG_LEVEL_ERROR,
  WARNING = SPDLOG_LEVEL_WARN,
  INFO = SPDLOG_LEVEL_INFO,
  DEBUG = SPDLOG_LEVEL_DEBUG,
};

void set_verbosity(Verbosity level);
void init(int argc, char** argv);

#define LOG_TRACE(...)                                                                       \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::trace, \
              __VA_ARGS__)
#ifndef NDEBUG
#define LOG_DEBUG(...)                                                                       \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::debug, \
              __VA_ARGS__)
#else
#define LOG_DEBUG(...) ((void)0)
#endif
#define LOG_INFO(...)                                                                       \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::info, \
              __VA_ARGS__)
#define LOG_WARN(...)                                                                       \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::warn, \
              __VA_ARGS__)
#define LOG_ERROR(...)                                                                     \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::err, \
              __VA_ARGS__)
#define LOG_FATAL(...)                                                                          \
  spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::critical, \
              __VA_ARGS__)
#define LOG_OFF(...) \
  do {               \
  } while (false)

}  // namespace oiseau::logging

template <typename T>
struct std::formatter<xt::xarray<T>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template <typename Context>
  auto format(xt::xarray<T> const& arr, Context& ctx) const {
    std::stringstream ss;
    ss << arr;
    return std::format_to(ctx.out(), "{}", ss.str());
  }
};
