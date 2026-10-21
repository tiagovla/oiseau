// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <ostream>

namespace oiseau::time {

/// Prints a bounded number of updates on one terminal line.
class ProgressBar {
 public:
  ProgressBar(std::ostream* output, std::size_t total) : m_output(output), m_total(total) {}

  void update(std::size_t completed) {
    if (m_output == nullptr || m_total == 0) return;

    constexpr std::size_t width = 40;
    const auto fraction = static_cast<long double>(completed) / static_cast<long double>(m_total);
    const auto filled = static_cast<std::size_t>(fraction * width);
    if (filled == m_last_filled && completed != m_total) return;
    m_last_filled = filled;

    *m_output << "\rIntegrating [";
    for (std::size_t column = 0; column < width; ++column) {
      *m_output << (column < filled ? '#' : ' ');
    }
    *m_output << "] " << static_cast<unsigned>(fraction * 100.0L) << "% (" << completed << '/'
              << m_total << ')';
    if (completed == m_total) *m_output << '\n';
    *m_output << std::flush;
  }

 private:
  std::ostream* m_output;
  std::size_t m_total;
  std::size_t m_last_filled = 41;
};

}  // namespace oiseau::time
