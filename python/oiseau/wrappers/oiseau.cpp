// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace oiseau::wrappers {
void mesh(py::module &m);
void io(py::module &m);
}  // namespace oiseau::wrappers

PYBIND11_MODULE(cpp, m) {
  m.doc() = "OISEAU Python interface";

  py::module mesh = m.def_submodule("mesh", "Mesh module");
  oiseau::wrappers::mesh(mesh);

  py::module io = m.def_submodule("io", "IO module");
  oiseau::wrappers::io(io);
}
