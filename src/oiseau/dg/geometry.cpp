// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/detail/cell_mapping.hpp"
#include "oiseau/dg/detail/geometry.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_element_type.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/utils/logging.hpp"
#include "xtensor/containers/xadapt.hpp"

namespace oiseau::dg {
namespace {

using detail::Matrix;
using detail::Point;

Matrix edge_matrix(const std::vector<Point>& points, const std::array<std::size_t, 4>& basis,
                   std::size_t dimension) {
  Matrix edges = xt::zeros<double>({dimension, dimension});
  for (std::size_t column = 0; column < dimension; ++column) {
    const auto edge = points[basis[column + 1]] - points[basis[0]];
    xt::view(edges, xt::all(), column) = xt::view(edge, xt::range(0, dimension));
  }
  return edges;
}

struct VertexBasis {
  std::array<std::size_t, 4> indices{};
  double determinant = 0.0;
};

void consider_basis_completions(const std::vector<Point>& points, std::size_t dimension,
                                std::size_t first, std::size_t second, VertexBasis& best) {
  const auto third_start = dimension > 2 ? second + 1 : 1;
  const auto third_end = dimension > 2 ? points.size() : 2;
  for (auto third = third_start; third < third_end; ++third) {
    const std::array<std::size_t, 4> indices = {0, first, second, third};
    const auto determinant = std::abs(xt::linalg::det(edge_matrix(points, indices, dimension)));
    if (determinant > best.determinant) {
      best = {.indices = indices, .determinant = determinant};
    }
  }
}

}  // namespace

namespace detail {

double triangle_area(const Point& a, const Point& b, const Point& c) {
  return xt::linalg::norm(xt::linalg::cross(b - a, c - a), 2) / 2.0;
}

std::array<double, 9> padded_matrix(const Matrix& matrix, std::size_t dimension) {
  std::array<double, 9> output{};
  auto padded = xt::adapt(output, std::array<std::size_t, 2>{3, 3});
  xt::view(padded, xt::range(0, dimension), xt::range(0, dimension)) = matrix;
  return output;
}

std::array<double, 3> point_array(const Point& point) {
  std::array<double, 3> output{};
  std::ranges::copy(std::views::take(point, 3), output.begin());
  return output;
}

bool is_degenerate(const Matrix& matrix, double determinant) {
  const auto scale = xt::amax(xt::abs(matrix))();
  if (!std::isfinite(determinant) || !std::isfinite(scale) || scale == 0.0) return true;
  const auto relative_threshold =
      100.0 * std::numeric_limits<double>::epsilon() * std::pow(scale, matrix.shape()[0]);
  return std::abs(determinant) <= relative_threshold;
}

double affine_tolerance(const std::vector<Point>& points, std::size_t dimension) {
  double coordinate_scale = 0.0;
  double cell_diameter = 0.0;

  for (auto point = points.cbegin(); point != points.cend(); ++point) {
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      coordinate_scale = std::max(coordinate_scale, std::abs((*point)(coordinate)));
    }
    cell_diameter =
        std::accumulate(std::next(point), points.cend(), cell_diameter,
                        [point](double diameter, const Point& other) {
                          return std::max(diameter, xt::linalg::norm(*point - other, 2));
                        });
  }

  constexpr double factor = 100.0 * std::numeric_limits<double>::epsilon();
  return factor * std::max(coordinate_scale, cell_diameter);
}

Point physical_vertex(const mesh::Geometry& geometry, std::size_t vertex, std::size_t dimension) {
  const auto coordinates = geometry.x();
  if (vertex >= coordinates.size() / geometry.dim()) {
    throw std::out_of_range("Cell references a vertex outside the mesh geometry");
  }

  Point point = xt::zeros<double>({dimension});
  for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
    point(coordinate) = coordinates[vertex * geometry.dim() + coordinate];
    if (!std::isfinite(point(coordinate))) {
      throw std::invalid_argument("Mesh coordinates must be finite");
    }
  }
  for (std::size_t coordinate = dimension; coordinate < geometry.dim(); ++coordinate) {
    const auto value = coordinates[vertex * geometry.dim() + coordinate];
    if (!std::isfinite(value)) throw std::invalid_argument("Mesh coordinates must be finite");
    if (value != 0.0) {
      throw std::invalid_argument("Geometry supports coordinates in the mesh dimensions only");
    }
  }
  return point;
}

std::array<std::size_t, 4> independent_vertices(const std::vector<Point>& reference_points,
                                                std::size_t dimension) {
  VertexBasis best;
  const auto second_end = dimension > 1 ? reference_points.size() : 2;
  for (std::size_t first = 1; first < reference_points.size(); ++first) {
    const auto second_start = dimension > 1 ? first + 1 : 1;
    for (auto second = second_start; second < second_end; ++second) {
      consider_basis_completions(reference_points, dimension, first, second, best);
    }
  }
  if (best.determinant <= std::numeric_limits<double>::epsilon()) {
    throw std::logic_error("Reference element has degenerate vertices");
  }
  return best.indices;
}

double face_measure(const std::vector<Point>& vertices, std::size_t face_dimension) {
  switch (face_dimension) {
  case 0:
    // A point face has measure 1 by convention.
    return 1.0;
  case 1:
    return xt::linalg::norm(vertices[1] - vertices[0], 2);
  case 2:
    if (vertices.size() == 3) {
      return triangle_area(vertices[0], vertices[1], vertices[2]);
    }
    return triangle_area(vertices[0], vertices[1], vertices[2]) +
           triangle_area(vertices[0], vertices[2], vertices[3]);
  default:
    throw std::logic_error("Unsupported face dimension");
  }
}

Point centroid(const std::vector<Point>& vertices) {
  Point output = xt::zeros<double>({vertices.front().size()});
  output = std::ranges::fold_left(vertices, std::move(output), [](Point sum, const Point& vertex) {
    sum += vertex;
    return sum;
  });
  output /= static_cast<double>(vertices.size());
  return output;
}

Point face_normal(const std::vector<Point>& vertices, const Point& cell_centroid,
                  std::size_t dimension) {
  Point normal = xt::zeros<double>({dimension});
  const auto face_centroid = centroid(vertices);
  if (dimension == 1) {
    normal = face_centroid - cell_centroid;
  } else if (dimension == 2) {
    const auto edge = vertices[1] - vertices[0];
    normal(0) = edge(1);
    normal(1) = -edge(0);
  } else {
    normal = xt::linalg::cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
  }

  const auto length = xt::linalg::norm(normal, 2);
  if (!std::isfinite(length) || length == 0.0) {
    throw std::invalid_argument("A cell face must have nonzero measure");
  }
  normal /= length;
  if (xt::linalg::vdot(normal, face_centroid - cell_centroid) < 0.0) {
    normal *= -1.0;
  }
  return normal;
}

}  // namespace detail

using namespace detail;

namespace {

struct ReferenceFaceGeometry {
  std::vector<std::size_t> vertices;
  double measure = 0.0;
};

struct ReferenceCellGeometry {
  // Points, basis indices, and face indices all use mesh-local vertex order.
  std::vector<Point> points;
  std::array<std::size_t, 4> basis_vertices{};
  Matrix inverse_edges;
  std::vector<ReferenceFaceGeometry> faces;
};

std::vector<ReferenceFaceGeometry> make_reference_faces(mesh::CellType cell_type,
                                                        const std::vector<Point>& points,
                                                        std::size_t face_dimension) {
  const auto local_faces = cell_type->get_entity_vertices(static_cast<int>(face_dimension));
  std::vector<ReferenceFaceGeometry> faces;
  faces.reserve(local_faces.size());
  for (const auto& local_face : local_faces) {
    std::vector<std::size_t> face_vertices;
    std::vector<Point> face_points;
    face_vertices.reserve(local_face.size());
    face_points.reserve(local_face.size());
    for (const auto vertex : local_face) {
      if (vertex < 0 || std::cmp_greater_equal(vertex, points.size())) {
        throw std::logic_error("Cell face references an invalid local vertex");
      }
      const auto index = static_cast<std::size_t>(vertex);
      face_vertices.push_back(index);
      face_points.push_back(points[index]);
    }
    const auto measure = face_measure(face_points, face_dimension);
    if (!std::isfinite(measure) || measure <= 0.0) {
      throw std::logic_error("Reference face has zero measure");
    }
    faces.push_back({.vertices = std::move(face_vertices), .measure = measure});
  }
  return faces;
}

ReferenceCellGeometry make_reference_geometry(mesh::CellType cell_type, std::size_t vertex_count,
                                              std::size_t dimension) {
  std::vector<std::size_t> local_indices(vertex_count);
  std::ranges::iota(local_indices, std::size_t{0});
  const auto reference_to_mesh = reference_ordered_vertices(cell_type, local_indices);
  const auto reference_element = nodal::get_ref_element(ref_type_for(cell_type->kind()), 1);
  if (reference_element->number_of_nodes() != vertex_count) {
    throw std::logic_error("Reference element and mesh cell vertex counts must match");
  }

  const auto& reference_nodes = reference_element->r();
  std::vector<Point> reference_points;
  reference_points.reserve(vertex_count);
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    Point point = xt::zeros<double>({dimension});
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      point(coordinate) =
          dimension == 1 ? reference_nodes(vertex) : reference_nodes(vertex, coordinate);
    }
    reference_points.push_back(std::move(point));
  }

  auto basis = independent_vertices(reference_points, dimension);
  const auto edges = edge_matrix(reference_points, basis, dimension);
  if (is_degenerate(edges, xt::linalg::det(edges))) {
    throw std::logic_error("Reference element has degenerate vertices");
  }
  for (std::size_t vertex = 0; vertex <= dimension; ++vertex) {
    basis[vertex] = reference_to_mesh[basis[vertex]];
  }

  ReferenceCellGeometry reference{
      .points = std::vector<Point>(vertex_count),
      .basis_vertices = basis,
      .inverse_edges = xt::linalg::inv(edges),
      .faces = {},
  };
  for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
    reference.points[reference_to_mesh[vertex]] = std::move(reference_points[vertex]);
  }

  reference.faces = make_reference_faces(cell_type, reference.points, dimension - 1);
  return reference;
}

std::size_t validated_mesh_dimension(const mesh::Topology& topology,
                                     const mesh::Geometry& geometry) {
  const auto dimension = static_cast<std::size_t>(topology.dimension());
  if (dimension < 1 || dimension > 3) {
    throw std::invalid_argument("Geometry supports meshes with dimensions 1 through 3");
  }
  if (geometry.dim() < dimension || geometry.x().size() % geometry.dim() != 0) {
    throw std::invalid_argument("Mesh coordinate data is inconsistent with its dimension");
  }
  if (topology.cell_types().size() != topology.conn().size()) {
    throw std::invalid_argument("Mesh cell types and connectivity counts must match");
  }
  return dimension;
}

Geometry::CellGeometry make_cell_geometry(const mesh::Geometry& geometry,
                                          const std::vector<std::size_t>& mesh_vertices,
                                          const ReferenceCellGeometry& reference,
                                          std::size_t dimension, std::size_t cell_index) {
  if (mesh_vertices.size() != reference.points.size()) {
    throw std::invalid_argument("Cell connectivity does not match its cell type");
  }

  std::vector<Point> physical_nodes;
  physical_nodes.reserve(mesh_vertices.size());
  std::ranges::transform(mesh_vertices, std::back_inserter(physical_nodes),
                         [&geometry, dimension](std::size_t vertex) {
                           return physical_vertex(geometry, vertex, dimension);
                         });
  const auto physical_edges = edge_matrix(physical_nodes, reference.basis_vertices, dimension);
  const auto jacobian = xt::linalg::dot(physical_edges, reference.inverse_edges);
  const auto determinant = xt::linalg::det(jacobian);
  if (is_degenerate(jacobian, determinant)) {
    throw std::invalid_argument("Cell vertices define a degenerate element");
  }
  const auto inverse_jacobian = xt::linalg::inv(jacobian);

  const auto tolerance = affine_tolerance(physical_nodes, dimension);
  const auto anchor = reference.basis_vertices[0];
  for (std::size_t vertex_index = 0; vertex_index < physical_nodes.size(); ++vertex_index) {
    const auto reference_offset = reference.points[vertex_index] - reference.points[anchor];
    const auto expected = xt::linalg::dot(jacobian, reference_offset);
    const auto actual = physical_nodes[vertex_index] - physical_nodes[anchor];
    if (xt::linalg::norm(expected - actual, 2) > tolerance) {
      throw std::invalid_argument("Geometry requires affine cell mappings");
    }
  }

  std::vector<Geometry::FaceGeometry> faces;
  faces.reserve(reference.faces.size());
  const auto cell_centroid = centroid(physical_nodes);
  for (const auto& reference_face : reference.faces) {
    std::vector<Point> physical_face;
    physical_face.reserve(reference_face.vertices.size());
    std::ranges::transform(
        reference_face.vertices, std::back_inserter(physical_face),
        [&physical_nodes](std::size_t vertex) -> const Point& { return physical_nodes[vertex]; });
    const auto measure = face_measure(physical_face, dimension - 1);
    if (!std::isfinite(measure) || measure <= 0.0) {
      throw std::invalid_argument("A cell face must have finite, nonzero measure");
    }
    faces.push_back(
        {.outward_normal = point_array(face_normal(physical_face, cell_centroid, dimension)),
         .measure = measure,
         .surface_jacobian = measure / reference_face.measure});
  }

  return {
      .mesh_cell_index = cell_index,
      .dimension = dimension,
      .jacobian = padded_matrix(jacobian, dimension),
      .inverse_jacobian = padded_matrix(inverse_jacobian, dimension),
      .determinant = determinant,
      .faces = std::move(faces),
  };
}

}  // namespace

Geometry::Geometry(const mesh::Mesh& mesh) : m_mesh(mesh) {
  const auto& topology = mesh.topology();
  const auto& geometry = mesh.geometry();
  const auto dimension = validated_mesh_dimension(topology, geometry);
  const auto cell_types = topology.cell_types();
  const auto connectivity = topology.conn();

  LOG_DEBUG("Computing DG geometry for a {}D mesh with {} cells", dimension, topology.n_cells());
  m_cells.resize(topology.n_cells());
  std::unordered_map<mesh::CellType, ReferenceCellGeometry> reference_geometry;
  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    const auto* cell_type = cell_types[cell_index];
    if (cell_type == nullptr) throw std::invalid_argument("Mesh cell type cannot be null");
    if (std::cmp_not_equal(cell_type->dimension(), dimension)) continue;

    const auto& mesh_vertices = connectivity[cell_index];
    auto [entry, inserted] = reference_geometry.try_emplace(cell_type);
    if (inserted) {
      entry->second = make_reference_geometry(cell_type, mesh_vertices.size(), dimension);
    }
    m_cells[cell_index] =
        make_cell_geometry(geometry, mesh_vertices, entry->second, dimension, cell_index);
  }
}

bool Geometry::has_cell(std::size_t mesh_cell_index) const {
  return mesh_cell_index < m_cells.size() && m_cells[mesh_cell_index].has_value();
}

const Geometry::CellGeometry& Geometry::cell(std::size_t mesh_cell_index) const {
  if (!has_cell(mesh_cell_index)) {
    throw std::out_of_range("No DG geometry is available for this mesh cell");
  }
  return *m_cells[mesh_cell_index];
}

}  // namespace oiseau::dg
