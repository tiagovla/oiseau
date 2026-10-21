// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>
#include <spdlog/details/log_msg.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>

#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "oiseau/utils/logging.hpp"

namespace {

class CollectingSink final : public spdlog::sinks::base_sink<std::mutex> {
 public:
  struct Record {
    spdlog::level::level_enum level;
    std::string message;
  };

  std::vector<Record> const& records() const { return records_; }

 protected:
  void sink_it_(spdlog::details::log_msg const& message) override {
    records_.push_back(
        {message.level, std::string(message.payload.data(), message.payload.size())});
  }

  void flush_() override {}

 private:
  std::vector<Record> records_;
};

class ScopedDefaultLogger {
 public:
  ScopedDefaultLogger() : original_(spdlog::default_logger()) {}

  ~ScopedDefaultLogger() { spdlog::set_default_logger(std::move(original_)); }

 private:
  std::shared_ptr<spdlog::logger> original_;
};

class ScopedLevel {
 public:
  ScopedLevel() : original_(spdlog::get_level()) {}

  ~ScopedLevel() { spdlog::set_level(original_); }

 private:
  spdlog::level::level_enum original_;
};

}  // namespace

TEST(test_utils_logging, set_verbosity_maps_levels) {
  ScopedLevel scoped_level;
  constexpr std::pair<oiseau::logging::Verbosity, spdlog::level::level_enum> levels[] = {
      {oiseau::logging::Verbosity::OFF, spdlog::level::off},
      {oiseau::logging::Verbosity::FATAL, spdlog::level::critical},
      {oiseau::logging::Verbosity::ERROR, spdlog::level::err},
      {oiseau::logging::Verbosity::WARNING, spdlog::level::warn},
      {oiseau::logging::Verbosity::INFO, spdlog::level::info},
      {oiseau::logging::Verbosity::DEBUG, spdlog::level::debug},
  };

  for (auto const& [verbosity, expected] : levels) {
    oiseau::logging::set_verbosity(verbosity);
    EXPECT_EQ(spdlog::get_level(), expected);
  }
}

TEST(test_utils_logging, init_sets_default_logger_and_info_level) {
  ScopedDefaultLogger scoped_default_logger;
  ScopedLevel scoped_level;
  char program_name[] = "oiseau_test_utils_logging";
  char* argv[] = {program_name, nullptr};

  oiseau::logging::init(1, argv);

  ASSERT_NE(spdlog::default_logger(), nullptr);
  EXPECT_EQ(spdlog::default_logger()->name(), "oiseau");
  EXPECT_EQ(spdlog::get("oiseau"), spdlog::default_logger());
  EXPECT_EQ(spdlog::get_level(), spdlog::level::info);
}

TEST(test_utils_logging, logging_macros_forward_levels_and_formatted_messages) {
  ScopedDefaultLogger scoped_default_logger;
  auto sink = std::make_shared<CollectingSink>();
  auto logger = std::make_shared<spdlog::logger>("oiseau_logging_test", sink);
  logger->set_level(spdlog::level::trace);
  spdlog::set_default_logger(logger);

  LOG_TRACE("trace {}", 1);
  LOG_DEBUG("debug {}", 2);
  LOG_INFO("info {}", 3);
  LOG_WARN("warning {}", 4);
  LOG_ERROR("error {}", 5);
  LOG_FATAL("fatal {}", 6);
  int off_argument_evaluations = 0;
  LOG_OFF("off {}", ++off_argument_evaluations);

  auto const& records = sink->records();
  ASSERT_EQ(records.size(), 6);
  EXPECT_EQ(records[0].level, spdlog::level::trace);
  EXPECT_EQ(records[0].message, "trace 1");
  EXPECT_EQ(records[1].level, spdlog::level::debug);
  EXPECT_EQ(records[1].message, "debug 2");
  EXPECT_EQ(records[2].level, spdlog::level::info);
  EXPECT_EQ(records[2].message, "info 3");
  EXPECT_EQ(records[3].level, spdlog::level::warn);
  EXPECT_EQ(records[3].message, "warning 4");
  EXPECT_EQ(records[4].level, spdlog::level::err);
  EXPECT_EQ(records[4].message, "error 5");
  EXPECT_EQ(records[5].level, spdlog::level::critical);
  EXPECT_EQ(records[5].message, "fatal 6");
  EXPECT_EQ(off_argument_evaluations, 0);
}
