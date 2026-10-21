// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "oiseau/dg/detail/cell_mapping.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/geometry.hpp"
#include "oiseau/mesh/topology.hpp"

namespace oiseau::dg {
namespace {

using Matrix = xt::xtensor<double, 2>;
using Point = xt::xtensor<double, 1>;

template <typename VectorExpression>
double vector_norm(const VectorExpression& vector) {
  return xt::linalg::norm(vector, 2);
}

Matrix zero_matrix(std::size_t dimension) { return xt::zeros<double>({dimension, dimension}); }

Point zero_point() { return xt::zeros<double>({3}); }

std::array<double, 9> padded_matrix(const Matrix& matrix, std::size_t dimension) {
  std::array<double, 9> output{};
  for (std::size_t row = 0; row < dimension; ++row) {
    for (std::size_t column = 0; column < dimension; ++column) {
      output[row * 3 + column] = matrix(row, column);
    }
  }
  return output;
}

std::array<double, 3> point_array(const Point& point) { return {point(0), point(1), point(2)}; }

nodal::RefElementType ref_type_for(mesh::CellKind kind) {
  switch (kind) {
  case mesh::CellKind::Interval:
    return nodal::RefElementType::Line;
  case mesh::CellKind::Triangle:
    return nodal::RefElementType::Triangle;
  case mesh::CellKind::Quadrilateral:
    return nodal::RefElementType::Quadrilateral;
  case mesh::CellKind::Tetrahedron:
    return nodal::RefElementType::Tetrahedron;
  case mesh::CellKind::Hexahedron:
    return nodal::RefElementType::Hexahedron;
  default:
    throw std::invalid_argument("DGGeometry does not support this mesh cell type");
  }
}

double matrix_scale(const Matrix& matrix, std::size_t dimension) {
  double scale = 0.0;
  for (std::size_t row = 0; row < dimension; ++row) {
    for (std::size_t column = 0; column < dimension; ++column) {
      scale = std::max(scale, std::abs(matrix(row, column)));
    }
  }
  return scale;
}

bool is_degenerate(const Matrix& matrix, double determinant, std::size_t dimension) {
  const auto scale = matrix_scale(matrix, dimension);
  if (!std::isfinite(determinant) || !std::isfinite(scale) || scale == 0.0) return true;
  const auto relative_threshold =
      100.0 * std::numeric_limits<double>::epsilon() * std::pow(scale, dimension);
  return std::abs(determinant) <= relative_threshold;
}

double affine_tolerance(const std::vector<Point>& points, std::size_t dimension) {
  double coordinate_scale = 0.0;
  for (const auto& point : points) {
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      coordinate_scale = std::max(coordinate_scale, std::abs(point(coordinate)));
    }
  }
  const auto max_distance_from = [&points](const Point& point) {
    return std::accumulate(points.begin(), points.end(), 0.0,
                           [&point](double distance, const Point& other) {
                             return std::max(distance, vector_norm(point - other));
                           });
  };
  const auto cell_diameter = std::accumulate(
      points.begin(), points.end(), 0.0, [&max_distance_from](double diameter, const Point& point) {
        return std::max(diameter, max_distance_from(point));
      });
  return 100.0 * std::numeric_limits<double>::epsilon() * std::max(coordinate_scale, cell_diameter);
}

Point physical_vertex(const mesh::Geometry& geometry, std::size_t vertex, std::size_t dimension) {
  const auto coordinates = geometry.x();
  if (vertex >= coordinates.size() / geometry.dim()) {
    throw std::out_of_range("Cell references a vertex outside the mesh geometry");
  }

  Point point = zero_point();
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
      throw std::invalid_argument("DGGeometry supports coordinates in the mesh dimensions only");
    }
  }
  return point;
}

double matrix_determinant(const Matrix& matrix) { return xt::linalg::det(matrix); }

Matrix inverse(const Matrix& matrix) { return xt::linalg::inv(matrix); }

std::array<std::size_t, 4> independent_vertices(const std::vector<Point>& reference_points,
                                                std::size_t dimension) {
  // Select independent reference edges from vertex zero. Some reference-node
  // orderings (notably tensor-product cells) do not put an independent basis
  // in their first dimension + 1 nodes.
  std::array<std::size_t, 4> result{};
  double largest_determinant = 0.0;
  for (std::size_t first = 1; first < reference_points.size(); ++first) {
    const auto second_end = dimension > 1 ? reference_points.size() : 2;
    const auto third_end = dimension > 2 ? reference_points.size() : 2;
    for (std::size_t second = 1; second < second_end; ++second) {
      for (std::size_t third = 1; third < third_end; ++third) {
        if ((dimension > 1 && second == first) ||
            (dimension > 2 && (third == first || third == second))) {
          continue;
        }
        const std::array<std::size_t, 3> indices = {first, second, third};
        auto edges = zero_matrix(dimension);
        for (std::size_t column = 0; column < dimension; ++column) {
          const auto edge = reference_points[indices[column]] - reference_points[0];
          for (std::size_t row = 0; row < dimension; ++row) {
            edges(row, column) = edge(row);
          }
        }
        const auto determinant = std::abs(matrix_determinant(edges));
        if (determinant > largest_determinant) {
          largest_determinant = determinant;
          result = {0, first, second, third};
        }
      }
    }
  }
  if (largest_determinant <= std::numeric_limits<double>::epsilon()) {
    throw std::logic_error("Reference element has degenerate vertices");
  }
  return result;
}

Point multiply(const Matrix& matrix, const Point& point, std::size_t dimension) {
  Point active_point = xt::zeros<double>({dimension});
  for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
    active_point(coordinate) = point(coordinate);
  }
  const auto product = xt::linalg::dot(matrix, active_point);
  Point output = zero_point();
  for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
    output(coordinate) = product(coordinate);
  }
  return output;
}

Matrix multiply(const Matrix& left, const Matrix& right) { return xt::linalg::dot(left, right); }

double face_measure(const std::vector<Point>& vertices, std::size_t dimension) {
  if (dimension == 1) return 1.0;
  if (dimension == 2) return vector_norm(vertices[1] - vertices[0]);
  if (vertices.size() == 3) {
    return vector_norm(xt::linalg::cross(vertices[1] - vertices[0], vertices[2] - vertices[0])) /
           2.0;
  }
  const auto first =
      vector_norm(xt::linalg::cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
  const auto second =
      vector_norm(xt::linalg::cross(vertices[2] - vertices[0], vertices[3] - vertices[0]));
  return (first + second) / 2.0;
}

Point centroid(const std::vector<Point>& vertices) {
  auto output = std::accumulate(vertices.begin(), vertices.end(), zero_point(),
                                [](Point sum, const Point& vertex) {
                                  sum += vertex;
                                  return sum;
                                });
  output /= static_cast<double>(vertices.size());
  return output;
}

Point face_normal(const std::vector<Point>& vertices, const Point& cell_centroid,
                  std::size_t dimension) {
  Point normal = zero_point();
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

  const auto length = vector_norm(normal);
  if (!std::isfinite(length) || length == 0.0) {
    throw std::invalid_argument("A cell face must have nonzero measure");
  }
  normal /= length;
  if (xt::linalg::vdot(normal, face_centroid - cell_centroid) < 0.0) {
    normal *= -1.0;
  }
  return normal;
}

}  // namespace

DGGeometry::DGGeometry(const mesh::Mesh& mesh) : m_mesh(mesh) {
  const auto& topology = mesh.topology();
  const auto& geometry = mesh.geometry();
  const auto cell_types = topology.cell_types();
  const auto connectivity = topology.conn();
  const auto dimension = static_cast<std::size_t>(topology.dimension());

  if (dimension < 1 || dimension > 3) {
    throw std::invalid_argument("DGGeometry supports meshes with dimensions 1 through 3");
  }
  if (geometry.dim() < dimension || geometry.x().size() % geometry.dim() != 0) {
    throw std::invalid_argument("Mesh coordinate data is inconsistent with its dimension");
  }
  if (cell_types.size() != connectivity.size()) {
    throw std::invalid_argument("Mesh cell types and connectivity counts must match");
  }

  m_cells.resize(topology.n_cells());
  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    const auto* cell_type = cell_types[cell_index];
    if (cell_type == nullptr) throw std::invalid_argument("Mesh cell type cannot be null");
    if (std::cmp_not_equal(cell_type->dimension(), dimension)) continue;

    const auto& mesh_vertices = connectivity[cell_index];
    auto vertices = detail::reference_ordered_vertices(cell_type, mesh_vertices);

    const auto reference_element = nodal::get_ref_element(ref_type_for(cell_type->kind()), 1);
    if (reference_element->number_of_nodes() != vertices.size()) {
      throw std::logic_error("Reference element and mesh cell vertex counts must match");
    }
    const auto& reference_nodes = reference_element->r();
    std::vector<Point> physical_nodes;
    std::vector<Point> reference_points;
    std::vector<Point> mesh_physical_nodes;
    std::vector<Point> mesh_reference_points(vertices.size());
    physical_nodes.reserve(vertices.size());
    reference_points.reserve(vertices.size());
    mesh_physical_nodes.reserve(vertices.size());
    std::ranges::transform(mesh_vertices, std::back_inserter(mesh_physical_nodes),
                           [&geometry, dimension](std::size_t vertex) {
                             return physical_vertex(geometry, vertex, dimension);
                           });
    for (std::size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {
      physical_nodes.push_back(physical_vertex(geometry, vertices[vertex_index], dimension));
      Point reference_point = zero_point();
      for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        reference_point(coordinate) = reference_nodes(vertex_index, coordinate);
      }
      reference_points.push_back(reference_point);
      const auto mesh_vertex = std::ranges::find(mesh_vertices, vertices[vertex_index]);
      mesh_reference_points[static_cast<std::size_t>(mesh_vertex - mesh_vertices.begin())] =
          reference_point;
    }

    auto reference_edges = zero_matrix(dimension);
    auto physical_edges = zero_matrix(dimension);
    const auto basis_vertices = independent_vertices(reference_points, dimension);
    for (std::size_t column = 0; column < dimension; ++column) {
      const auto reference_edge =
          reference_points[basis_vertices[column + 1]] - reference_points[basis_vertices[0]];
      const auto physical_edge =
          physical_nodes[basis_vertices[column + 1]] - physical_nodes[basis_vertices[0]];
      for (std::size_t row = 0; row < dimension; ++row) {
        reference_edges(row, column) = reference_edge(row);
        physical_edges(row, column) = physical_edge(row);
      }
    }

    const auto reference_determinant = matrix_determinant(reference_edges);
    if (is_degenerate(reference_edges, reference_determinant, dimension)) {
      throw std::logic_error("Reference element has degenerate vertices");
    }
    const auto inverse_reference = inverse(reference_edges);
    const auto jacobian = multiply(physical_edges, inverse_reference);
    const auto determinant = matrix_determinant(jacobian);
    if (is_degenerate(jacobian, determinant, dimension)) {
      throw std::invalid_argument("Cell vertices define a degenerate element");
    }
    const auto inverse_jacobian = inverse(jacobian);

    const auto tolerance = affine_tolerance(physical_nodes, dimension);
    for (std::size_t vertex_index = 0; vertex_index < physical_nodes.size(); ++vertex_index) {
      const auto reference_offset =
          reference_points[vertex_index] - reference_points[basis_vertices[0]];
      const auto expected = multiply(jacobian, reference_offset, dimension);
      const auto actual = physical_nodes[vertex_index] - physical_nodes[basis_vertices[0]];
      if (vector_norm(expected - actual) > tolerance) {
        throw std::invalid_argument("DGGeometry requires affine cell mappings");
      }
    }

    const auto local_faces = cell_type->get_entity_vertices(static_cast<int>(dimension - 1));
    std::vector<FaceGeometry> faces;
    faces.reserve(local_faces.size());
    const auto cell_centroid = centroid(mesh_physical_nodes);
    for (const auto& local_face : local_faces) {
      std::vector<Point> physical_face;
      std::vector<Point> reference_face;
      physical_face.reserve(local_face.size());
      reference_face.reserve(local_face.size());
      for (const auto local_vertex : local_face) {
        if (local_vertex < 0 ||
            static_cast<std::size_t>(local_vertex) >= mesh_physical_nodes.size()) {
          throw std::logic_error("Cell face references an invalid local vertex");
        }
        physical_face.push_back(mesh_physical_nodes[static_cast<std::size_t>(local_vertex)]);
        reference_face.push_back(mesh_reference_points[static_cast<std::size_t>(local_vertex)]);
      }
      const auto measure = face_measure(physical_face, dimension);
      const auto reference_measure = face_measure(reference_face, dimension);
      if (!std::isfinite(measure) || measure <= 0.0) {
        throw std::invalid_argument("A cell face must have finite, nonzero measure");
      }
      if (!std::isfinite(reference_measure) || reference_measure <= 0.0) {
        throw std::logic_error("Reference face has zero measure");
      }
      faces.push_back({point_array(face_normal(physical_face, cell_centroid, dimension)), measure,
                       measure / reference_measure});
    }

    m_cells[cell_index] = CellGeometry{
        .mesh_cell_index = cell_index,
        .dimension = dimension,
        .jacobian = padded_matrix(jacobian, dimension),
        .inverse_jacobian = padded_matrix(inverse_jacobian, dimension),
        .determinant = determinant,
        .faces = std::move(faces),
    };
  }
}

bool DGGeometry::has_cell(std::size_t mesh_cell_index) const {
  return mesh_cell_index < m_cells.size() && m_cells[mesh_cell_index].has_value();
}

const DGGeometry::CellGeometry& DGGeometry::cell(std::size_t mesh_cell_index) const {
  if (!has_cell(mesh_cell_index)) {
    throw std::out_of_range("No DG geometry is available for this mesh cell");
  }
  return *m_cells[mesh_cell_index];
}

}  // namespace oiseau::dg
