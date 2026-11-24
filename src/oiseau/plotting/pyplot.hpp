// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <memory>
#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xexpression.hpp>
#include <xtensor/core/xtensor_forward.hpp>

namespace py = nanobind;
using namespace py::literals;

namespace nanobind::detail {

template <typename E>
  requires xt::is_xexpression<E>::value
struct type_caster<E> {
  using value_type = typename E::value_type;
  using array_type = xt::xarray<value_type>;

  NB_TYPE_CASTER(array_type, const_name("numpy.ndarray"));
  static handle from_cpp(const E& expr, rv_policy, cleanup_list*) {
    auto data = std::make_unique<array_type>(expr);
    py::capsule owner(data.get(),
                      [](void* pointer) noexcept { delete static_cast<array_type*>(pointer); });
    auto* array = data.release();
    return py::ndarray<py::numpy, value_type>(array->data(), array->dimension(),
                                              array->shape().data(), owner)
        .cast()
        .release();
  }
};

}  // namespace nanobind::detail

namespace plt {

class scoped_interpreter {
 public:
  NB_INLINE scoped_interpreter() {
    // Register nanobind state in the calling binary as well as the plotting library.
    if (!py::register_module(initialize())) fail_initialization();
  }
  ~scoped_interpreter();

  scoped_interpreter(const scoped_interpreter&) = delete;
  scoped_interpreter& operator=(const scoped_interpreter&) = delete;
  scoped_interpreter(scoped_interpreter&& other) noexcept
      : m_active(std::exchange(other.m_active, false)) {}

 private:
  static NB_EXPORT py::module_ initialize();
  [[noreturn]] static void fail_initialization();
  bool m_active = true;
};

inline auto plt() { return py::module_::import_("matplotlib.pyplot"); }

#define DEFINE_PYPLOT_FUNC(name)                           \
  template <typename... Args>                              \
  auto name(Args&&... args) {                              \
    return plt().attr(#name)(std::forward<Args>(args)...); \
  }

#define DEFINE_PYPLOT_FUNC_WRAPPED(WrapperClass, name)                   \
  template <typename... Args>                                            \
  WrapperClass name(Args&&... args) {                                    \
    return WrapperClass(plt().attr(#name)(std::forward<Args>(args)...)); \
  }

#define DEFINE_PY_CLASS_METHOD_AUTO(name)                  \
  template <typename... Args>                              \
  auto name(Args&&... args) {                              \
    return m_obj.attr(#name)(std::forward<Args>(args)...); \
  }

#define DEFINE_PY_CLASS_METHOD_WRAPPED(WrapperClass, name)               \
  template <typename... Args>                                            \
  WrapperClass name(Args&&... args) {                                    \
    return WrapperClass(m_obj.attr(#name)(std::forward<Args>(args)...)); \
  }

// --- Classes using the simplified macros ---

class NB_EXPORT AxesSubPlot {
 public:
  explicit AxesSubPlot(py::object axes) : m_obj(std::move(axes)) {}

  DEFINE_PY_CLASS_METHOD_AUTO(scatter)
  DEFINE_PY_CLASS_METHOD_AUTO(triplot)
  DEFINE_PY_CLASS_METHOD_AUTO(tripcolor)
  DEFINE_PY_CLASS_METHOD_AUTO(plot)
  DEFINE_PY_CLASS_METHOD_AUTO(quiver)
  DEFINE_PY_CLASS_METHOD_AUTO(text)
  DEFINE_PY_CLASS_METHOD_AUTO(set_xlabel)
  DEFINE_PY_CLASS_METHOD_AUTO(set_ylabel)
  DEFINE_PY_CLASS_METHOD_AUTO(set_zlabel)
  DEFINE_PY_CLASS_METHOD_AUTO(set_title)
  DEFINE_PY_CLASS_METHOD_AUTO(set_aspect)
  DEFINE_PY_CLASS_METHOD_AUTO(grid)
  DEFINE_PY_CLASS_METHOD_AUTO(legend)
  DEFINE_PY_CLASS_METHOD_AUTO(set_xlim)
  DEFINE_PY_CLASS_METHOD_AUTO(set_ylim)
  DEFINE_PY_CLASS_METHOD_AUTO(set_zlim)

 private:
  py::object m_obj;
};

class NB_EXPORT Figure {
 public:
  explicit Figure(py::object fig) : m_obj(std::move(fig)) {}

  DEFINE_PY_CLASS_METHOD_WRAPPED(AxesSubPlot, add_subplot)
  DEFINE_PY_CLASS_METHOD_AUTO(savefig)
  DEFINE_PY_CLASS_METHOD_AUTO(show)

 private:
  py::object m_obj;
};

// --- Free functions using the simplified macros ---

DEFINE_PYPLOT_FUNC(scatter)
DEFINE_PYPLOT_FUNC(show)
DEFINE_PYPLOT_FUNC(plot)
DEFINE_PYPLOT_FUNC(xlabel)
DEFINE_PYPLOT_FUNC(ylabel)
DEFINE_PYPLOT_FUNC(title)
DEFINE_PYPLOT_FUNC(grid)
DEFINE_PYPLOT_FUNC(legend)
DEFINE_PYPLOT_FUNC(savefig)

DEFINE_PYPLOT_FUNC_WRAPPED(Figure, figure)

template <typename... Args>
std::pair<Figure, AxesSubPlot> subplots(Args&&... args) {
  auto res = py::cast<py::tuple>(plt().attr("subplots")(std::forward<Args>(args)...));
  return {Figure(res[0]), AxesSubPlot(res[1])};
}

// --- Undefine all macros ---

#undef DEFINE_PYPLOT_FUNC
#undef DEFINE_PYPLOT_FUNC_WRAPPED
#undef DEFINE_PY_CLASS_METHOD_AUTO
#undef DEFINE_PY_CLASS_METHOD_WRAPPED

}  // namespace plt
