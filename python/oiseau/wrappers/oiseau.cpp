// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <nanobind/nanobind.h>

namespace py = nanobind;

namespace oiseau::wrappers {
void mesh(py::module_& m);
void io(py::module_& m);
void dg(py::module_& m);
}  // namespace oiseau::wrappers

NB_MODULE(cpp, m) {
  m.attr("__doc__") = "OISEAU Python interface";

  py::module_ mesh = m.def_submodule("mesh", "Mesh module");
  oiseau::wrappers::mesh(mesh);

  py::module_ io = m.def_submodule("io", "IO module");
  oiseau::wrappers::io(io);

  py::module_ dg = m.def_submodule("dg", "Discontinuous Galerkin solvers");
  oiseau::wrappers::dg(dg);
}
