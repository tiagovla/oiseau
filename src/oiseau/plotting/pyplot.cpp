// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/plotting/pyplot.hpp"

#include <Python.h>
#include <cpython/initconfig.h>
#include <import.h>
#include <nanobind/nanobind.h>
#include <pylifecycle.h>
#include <pythonrun.h>

#include <stdexcept>

NB_MODULE(_oiseau_plotting, module) { (void)module; }

namespace plt {

py::module_ scoped_interpreter::initialize() {
  if (Py_IsInitialized()) throw std::runtime_error("The Python interpreter is already running");

  // Importing a built-in module initializes nanobind in the embedded interpreter.
  static const int registered = PyImport_AppendInittab("_oiseau_plotting", PyInit__oiseau_plotting);
  if (registered == -1) throw std::runtime_error("Could not register the plotting module");

  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  config.parse_argv = 0;
  const auto status = Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);
  if (PyStatus_Exception(status)) {
    throw std::runtime_error(status.err_msg ? status.err_msg : "Could not initialize Python");
  }

  auto module = py::steal<py::module_>(PyImport_ImportModule("_oiseau_plotting"));
  if (!module.is_valid()) fail_initialization();
  return module;
}

void scoped_interpreter::fail_initialization() {
  PyErr_Print();
  Py_FinalizeEx();
  throw std::runtime_error("Could not initialize the plotting module");
}

scoped_interpreter::~scoped_interpreter() {
  if (m_active) Py_FinalizeEx();
}

}  // namespace plt
