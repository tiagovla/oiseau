// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstddef>
#include <iostream>
#include <span>
#include <vector>

#include "oiseau/dg/nodal/element.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace oiseau::dg {

class DGSpace {
 public:
  DGSpace(const DGSpace& V) = delete;
  DGSpace(const mesh::Mesh& mesh, const std::vector<unsigned>& orders);
  DGSpace(DGSpace&& V) = default;
  virtual ~DGSpace() = default;
  DGSpace& operator=(const DGSpace& V) = delete;
  DGSpace& operator=(DGSpace&& V) = delete;

  inline const mesh::Mesh& mesh() const { return m_mesh; };
  std::span<const nodal::Element> elements() const;
  std::span<const unsigned> orders() const;
  std::span<const std::size_t> cell_indices() const;

 private:
  const mesh::Mesh& m_mesh;
  std::vector<nodal::Element> m_elements;
  const std::vector<unsigned> m_orders;
  std::vector<std::size_t> m_cell_indices;
};

}  // namespace oiseau::dg
