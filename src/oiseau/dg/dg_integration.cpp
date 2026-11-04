// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_integration.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/views/xslice.hpp>

#include "oiseau/dg/detail/cell_mapping.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/utils/integration.hpp"

namespace oiseau::dg {
namespace {

using Matrix = DGIntegration::Matrix;
using Point2 = std::array<double, 2>;

struct ReferenceOperators {
  std::shared_ptr<nodal::RefElement> element;
  Matrix mass;
  Matrix inverse_mass;
  Matrix stiffness_r;
  Matrix stiffness_s;
  Matrix inverse_vandermonde;
};

nodal::RefElementType reference_type(mesh::CellKind kind) {
  switch (kind) {
  case mesh::CellKind::Triangle:
    return nodal::RefElementType::Triangle;
  case mesh::CellKind::Quadrilateral:
    return nodal::RefElementType::Quadrilateral;
  default:
    throw std::invalid_argument("DGIntegration supports triangles and quadrilaterals only");
  }
}

ReferenceOperators make_reference_operators(mesh::CellKind kind, unsigned order) {
  auto element = nodal::get_ref_element(reference_type(kind), order);
  const auto vandermonde = element->v();
  const auto inverse_vandermonde = xt::linalg::inv(vandermonde);

  Matrix mass = xt::linalg::dot(xt::transpose(inverse_vandermonde), inverse_vandermonde);
  Matrix inverse_mass = xt::linalg::dot(vandermonde, xt::transpose(vandermonde));

  const auto derivatives = element->d();
  const Matrix derivative_r = xt::view(derivatives, xt::all(), xt::all(), 0);
  const Matrix derivative_s = xt::view(derivatives, xt::all(), xt::all(), 1);
  Matrix stiffness_r = xt::linalg::dot(xt::transpose(derivative_r), mass);
  Matrix stiffness_s = xt::linalg::dot(xt::transpose(derivative_s), mass);

  return {.element = std::move(element),
          .mass = std::move(mass),
          .inverse_mass = std::move(inverse_mass),
          .stiffness_r = std::move(stiffness_r),
          .stiffness_s = std::move(stiffness_s),
          .inverse_vandermonde = std::move(inverse_vandermonde)};
}

std::vector<Point2> reference_vertex_coordinates(mesh::CellType type,
                                                 const std::vector<std::size_t>& mesh_vertices) {
  const auto vertices = detail::reference_ordered_vertices(type, mesh_vertices);
  const auto reference = nodal::get_ref_element(reference_type(type->kind()), 1);
  std::vector<Point2> output(mesh_vertices.size());
  for (std::size_t reference_vertex = 0; reference_vertex < vertices.size(); ++reference_vertex) {
    const auto position = std::ranges::find(mesh_vertices, vertices[reference_vertex]);
    if (position == mesh_vertices.end()) {
      throw std::logic_error("Reference vertex is missing from mesh connectivity");
    }
    const auto mesh_vertex = static_cast<std::size_t>(position - mesh_vertices.begin());
    output[mesh_vertex] = {reference->r()(reference_vertex, 0),
                           reference->r()(reference_vertex, 1)};
  }
  return output;
}

Matrix edge_points(const Point2& first, const Point2& second, const Matrix& line_nodes) {
  Matrix points = Matrix::from_shape({line_nodes.size(), 2});
  for (std::size_t q = 0; q < line_nodes.size(); ++q) {
    const auto coordinate = line_nodes(q);
    const auto first_weight = 0.5 * (1.0 - coordinate);
    const auto second_weight = 0.5 * (1.0 + coordinate);
    for (std::size_t d = 0; d < 2; ++d) {
      points(q, d) = first_weight * first[d] + second_weight * second[d];
    }
  }
  return points;
}

Matrix interpolation_on_points(const ReferenceOperators& operators, const Matrix& points) {
  const auto vandermonde = operators.element->vandermonde(points);
  return xt::linalg::dot(vandermonde, operators.inverse_vandermonde);
}

}  // namespace

DGIntegration::DGIntegration(const mesh::Mesh& mesh, const DGSpace& space,
                             const DGConnectivity& connectivity, const DGGeometry& geometry) {
  if (&space.mesh() != &mesh || &connectivity.mesh() != &mesh || &geometry.mesh() != &mesh) {
    throw std::invalid_argument("DGIntegration inputs must refer to the same mesh");
  }

  const auto groups = space.element_groups();
  const auto cell_types = mesh.topology().cell_types();
  const auto connectivity_data = mesh.topology().conn();
  const auto cell_count = mesh.topology().n_cells();

  std::vector<ReferenceOperators> reference_operators;
  std::vector<std::size_t> cell_group(cell_count, DGConnectivity::no_neighbor);
  std::vector<std::size_t> cell_group_row(cell_count, 0);
  reference_operators.reserve(groups.size());
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    const auto& group = groups[group_index];
    reference_operators.push_back(make_reference_operators(group.kind, group.order));
    for (std::size_t row = 0; row < group.cell_indices.size(); ++row) {
      const auto cell_index = group.cell_indices[row];
      cell_group[cell_index] = group_index;
      cell_group_row[cell_index] = row;
    }
  }

  m_cells.resize(cell_count);
  for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
    if (cell_group[cell_index] == DGConnectivity::no_neighbor) continue;

    const auto group_index = cell_group[cell_index];
    const auto row = cell_group_row[cell_index];
    const auto& reference = reference_operators[group_index];
    const auto& cell_geometry = geometry.cell(cell_index);
    const auto determinant = std::abs(cell_geometry.determinant);
    auto& cell_rule = m_cells[cell_index];
    cell_rule.group_index = group_index;
    cell_rule.group_cell_index = row;
    cell_rule.mass = determinant * reference.mass;
    cell_rule.inverse_mass = reference.inverse_mass / determinant;
    cell_rule.stiffness_x =
        determinant * (cell_geometry.inverse_jacobian[0] * reference.stiffness_r +
                       cell_geometry.inverse_jacobian[3] * reference.stiffness_s);
    cell_rule.stiffness_y =
        determinant * (cell_geometry.inverse_jacobian[1] * reference.stiffness_r +
                       cell_geometry.inverse_jacobian[4] * reference.stiffness_s);

    const auto* cell_type = cell_types[cell_index];
    const auto local_faces = cell_type->get_entity_vertices(1);
    const auto face_connections = connectivity.faces(cell_index);
    cell_rule.faces.reserve(local_faces.size());
    for (std::size_t face_index = 0; face_index < local_faces.size(); ++face_index) {
      const auto& connection = face_connections[face_index];
      auto& face_rule = cell_rule.faces.emplace_back();
      face_rule.outward_normal = {cell_geometry.faces[face_index].outward_normal[0],
                                  cell_geometry.faces[face_index].outward_normal[1]};
      face_rule.neighbor_cell = connection.neighbor_cell;

      auto neighbor_group = group_index;
      auto quadrature_order = groups[group_index].order;
      if (!connection.is_boundary()) {
        neighbor_group = cell_group.at(connection.neighbor_cell);
        if (neighbor_group == DGConnectivity::no_neighbor) {
          throw std::invalid_argument("DGIntegration neighbor cell is absent from DGSpace");
        }
        quadrature_order = std::max(quadrature_order, groups[neighbor_group].order);
      }
      if (quadrature_order > (static_cast<unsigned>(std::numeric_limits<int>::max()) - 1U) / 2U) {
        throw std::invalid_argument("Polynomial order is too large for DG face quadrature");
      }
      auto [line_nodes, line_weights] =
          utils::integration::quadrature(static_cast<int>(2U * quadrature_order + 1U));

      const auto& mesh_vertices = connectivity_data[cell_index];
      const auto reference_vertices = reference_vertex_coordinates(cell_type, mesh_vertices);
      const auto first_local_vertex = static_cast<std::size_t>(local_faces[face_index][0]);
      const auto second_local_vertex = static_cast<std::size_t>(local_faces[face_index][1]);
      const auto local_points = edge_points(reference_vertices[first_local_vertex],
                                            reference_vertices[second_local_vertex], line_nodes);
      face_rule.minus_trace = interpolation_on_points(reference, local_points);

      const auto face_measure = cell_geometry.faces[face_index].measure;
      face_rule.physical_weights.resize(line_weights.size());
      for (std::size_t q = 0; q < line_weights.size(); ++q) {
        face_rule.physical_weights[q] = 0.5 * face_measure * line_weights(q);
      }

      if (connection.is_boundary()) continue;

      const auto neighbor_index = connection.neighbor_cell;
      const auto* neighbor_type = cell_types[neighbor_index];
      const auto& neighbor_mesh_vertices = connectivity_data[neighbor_index];
      const auto neighbor_reference_vertices =
          reference_vertex_coordinates(neighbor_type, neighbor_mesh_vertices);
      const auto neighbor_local_faces = neighbor_type->get_entity_vertices(1);
      const auto& neighbor_face_vertices = neighbor_local_faces[connection.neighbor_face];
      if (connection.vertex_permutation.size() != 2 || neighbor_face_vertices.size() != 2) {
        throw std::logic_error("DGIntegration expects two vertices on every face");
      }
      const auto neighbor_first =
          static_cast<std::size_t>(neighbor_face_vertices[connection.vertex_permutation[0]]);
      const auto neighbor_second =
          static_cast<std::size_t>(neighbor_face_vertices[connection.vertex_permutation[1]]);
      const auto neighbor_points =
          edge_points(neighbor_reference_vertices[neighbor_first],
                      neighbor_reference_vertices[neighbor_second], line_nodes);
      face_rule.plus_trace =
          interpolation_on_points(reference_operators[neighbor_group], neighbor_points);
    }
  }
}

const DGIntegration::CellRule& DGIntegration::cell(std::size_t mesh_cell_index) const {
  return m_cells.at(mesh_cell_index);
}

}  // namespace oiseau::dg
