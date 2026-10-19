// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/nodal/element.hpp"

#include <memory>
#include <utility>
#include <xtensor/core/xtensor_forward.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"

namespace oiseau::dg::nodal {

Element::Element(std::shared_ptr<RefElement> ref_elem, xt::xarray<double> nodes)
    : m_reference(std::move(ref_elem)), m_nodes(std::move(nodes)) {}
const RefElement& Element::reference() const { return *m_reference; }
unsigned Element::order() const { return this->m_reference->order(); }
const xt::xarray<double>& Element::nodes() const { return this->m_nodes; }
}  // namespace oiseau::dg::nodal
