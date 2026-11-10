// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"

namespace oiseau::dg::detail {

using Matrix = xt::xtensor<double, 2>;
using Point = xt::xtensor<double, 1>;

double vector_norm(const Point& vector);
std::array<double, 9> padded_matrix(const Matrix& matrix, std::size_t dimension);
std::array<double, 3> point_array(const Point& point);
nodal::RefElementType ref_type_for(mesh::CellKind kind);
double matrix_scale(const Matrix& matrix, std::size_t dimension);
bool is_degenerate(const Matrix& matrix, double determinant, std::size_t dimension);
double affine_tolerance(const std::vector<Point>& points, std::size_t dimension);
Point physical_vertex(const mesh::Geometry& geometry, std::size_t vertex, std::size_t dimension);
std::array<std::size_t, 4> independent_vertices(const std::vector<Point>& reference_points,
                                                std::size_t dimension);
double face_measure(const std::vector<Point>& vertices, std::size_t dimension);
Point centroid(const std::vector<Point>& vertices);
Point face_normal(const std::vector<Point>& vertices, const Point& cell_centroid,
                  std::size_t dimension);

}  // namespace oiseau::dg::detail
