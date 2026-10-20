// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace oiseau::dg {

/// A named scalar unknown used to declare a differential equation system.
///
/// Copies of a Field keep the same identity, so equation expressions remain
/// connected to the declared unknown even when handles are copied.
class Field {
 public:
  explicit Field(std::string name) : m_identity(std::make_shared<const Identity>(std::move(name))) {
    if (m_identity->name.empty()) {
      throw std::invalid_argument("Field name cannot be empty");
    }
  }

  const std::string& name() const { return m_identity->name; }

  bool same_as(const Field& other) const { return m_identity == other.m_identity; }

 private:
  struct Identity {
    explicit Identity(std::string field_name) : name(std::move(field_name)) {}
    std::string name;
  };

  std::shared_ptr<const Identity> m_identity;
};

}  // namespace oiseau::dg
