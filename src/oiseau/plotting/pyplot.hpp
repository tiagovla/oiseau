// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/detail/descr.h>
#include <pybind11/embed.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <utility>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xexpression.hpp>
#include <xtensor/core/xtensor_forward.hpp>

namespace py = pybind11;
using namespace py::literals;

namespace pybind11::detail {

template <typename T>
struct type_caster<xt::xarray<T>> {
 public:
  PYBIND11_TYPE_CASTER(xt::xarray<T>, const_name("xarray"));
  static handle cast(xt::xarray<T>& src, return_value_policy, handle) {
    auto v = py::array(src.size(), src.data(), py::capsule([]() {}));
    v.resize(src.shape());
    return v.release();
  }
};

template <typename E>
  requires xt::is_xexpression<E>::value
struct type_caster<E> {
  using value_type = typename E::value_type;
  using array_type = xt::xarray<value_type>;

  PYBIND11_TYPE_CASTER(array_type, _("xexpression"));
  static handle cast(const E& expr, return_value_policy, handle) {
    xt::xarray<value_type> tmp = xt::eval(expr);
    auto* data_holder = new xt::xarray<double>(std::move(tmp));
    py::capsule owner_capsule(data_holder,
                              [](void* p) { delete static_cast<xt::xarray<value_type>*>(p); });
    py::array result(data_holder->shape(), data_holder->data(), owner_capsule);
    return result.release();
  }
};

}  // namespace pybind11::detail

namespace plt {

using scoped_interpreter = py::scoped_interpreter;

inline auto plt() { return py::module_::import("matplotlib.pyplot"); }

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

class PYBIND11_EXPORT AxesSubPlot {
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

class PYBIND11_EXPORT Figure {
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
  py::object result = plt().attr("subplots")(std::forward<Args>(args)...);
  auto res = result.cast<py::tuple>();
  return {Figure(res[0]), AxesSubPlot(res[1])};
}

// --- Undefine all macros ---

#undef DEFINE_PYPLOT_FUNC
#undef DEFINE_PYPLOT_FUNC_WRAPPED
#undef DEFINE_PY_CLASS_METHOD_AUTO
#undef DEFINE_PY_CLASS_METHOD_WRAPPED

}  // namespace plt
