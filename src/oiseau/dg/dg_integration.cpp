// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_integration.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <tuple>
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
  std::shared_ptr<const DGIntegration::ReferenceRule> volume;
  Matrix inverse_vandermonde;
  std::vector<Point2> vertices;
  std::vector<std::vector<int>> faces;
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

  Matrix weak_derivative_r = xt::linalg::dot(inverse_mass, stiffness_r);
  Matrix weak_derivative_s = xt::linalg::dot(inverse_mass, stiffness_s);
  auto volume = std::make_shared<const DGIntegration::ReferenceRule>(
      DGIntegration::ReferenceRule{.mass = std::move(mass),
                                   .inverse_mass = std::move(inverse_mass),
                                   .stiffness_r = std::move(stiffness_r),
                                   .stiffness_s = std::move(stiffness_s),
                                   .weak_derivative_r = std::move(weak_derivative_r),
                                   .weak_derivative_s = std::move(weak_derivative_s)});
  return {.element = std::move(element),
          .volume = std::move(volume),
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
    auto& reference =
        reference_operators.emplace_back(make_reference_operators(group.kind, group.order));
    const auto first_cell = group.cell_indices.front();
    reference.vertices =
        reference_vertex_coordinates(cell_types[first_cell], connectivity_data[first_cell]);
    reference.faces = cell_types[first_cell]->get_entity_vertices(1);
    for (std::size_t row = 0; row < group.cell_indices.size(); ++row) {
      const auto cell_index = group.cell_indices[row];
      cell_group[cell_index] = group_index;
      cell_group_row[cell_index] = row;
    }
  }

  std::map<unsigned, std::shared_ptr<const QuadratureRule>> quadratures;
  using FaceKey = std::tuple<std::size_t, std::size_t, unsigned, bool>;
  std::map<FaceKey, std::shared_ptr<const ReferenceFaceRule>> face_rules;
  const auto face_reference = [&](std::size_t group_index, std::size_t face_index, unsigned order,
                                  bool reversed) {
    const FaceKey key{group_index, face_index, order, reversed};
    if (const auto found = face_rules.find(key); found != face_rules.end()) return found->second;

    auto& quadrature = quadratures[order];
    if (!quadrature) {
      if (order > (static_cast<unsigned>(std::numeric_limits<int>::max()) - 1U) / 2U) {
        throw std::invalid_argument("Polynomial order is too large for DG face quadrature");
      }
      auto [nodes, weights] = utils::integration::quadrature(static_cast<int>(2U * order + 1U));
      quadrature = std::make_shared<const QuadratureRule>(
          QuadratureRule{.nodes = std::move(nodes), .weights = std::move(weights)});
    }
    const auto& reference = reference_operators[group_index];
    const auto& vertices = reference.faces.at(face_index);
    if (vertices.size() != 2) {
      throw std::logic_error("DGIntegration expects two vertices on every face");
    }
    const auto first = static_cast<std::size_t>(vertices[reversed ? 1 : 0]);
    const auto second = static_cast<std::size_t>(vertices[reversed ? 0 : 1]);
    const auto points =
        edge_points(reference.vertices[first], reference.vertices[second], quadrature->nodes);
    Matrix trace = interpolation_on_points(reference, points);
    Matrix lift = xt::linalg::dot(trace, xt::transpose(reference.volume->inverse_mass));
    for (std::size_t q = 0; q < quadrature->weights.size(); ++q) {
      for (std::size_t node = 0; node < lift.shape()[1]; ++node) {
        lift(q, node) *= quadrature->weights(q);
      }
    }
    auto rule = std::make_shared<const ReferenceFaceRule>(ReferenceFaceRule{
        .quadrature = quadrature, .trace = std::move(trace), .lift = std::move(lift)});
    face_rules.emplace(key, rule);
    return rule;
  };

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
    cell_rule.reference = reference.volume;
    cell_rule.volume_jacobian = determinant;
    cell_rule.inverse_jacobian = {
        cell_geometry.inverse_jacobian[0], cell_geometry.inverse_jacobian[1],
        cell_geometry.inverse_jacobian[3], cell_geometry.inverse_jacobian[4]};

    const auto nodes = cell_rule.node_count();
    cell_rule.weak_derivative_x.resize(nodes * nodes);
    cell_rule.weak_derivative_y.resize(nodes * nodes);
    for (std::size_t test = 0; test < nodes; ++test) {
      for (std::size_t trial = 0; trial < nodes; ++trial) {
        const auto r = reference.volume->weak_derivative_r(test, trial);
        const auto s = reference.volume->weak_derivative_s(test, trial);
        cell_rule.weak_derivative_x[test * nodes + trial] =
            cell_rule.inverse_jacobian[0] * r + cell_rule.inverse_jacobian[2] * s;
        cell_rule.weak_derivative_y[test * nodes + trial] =
            cell_rule.inverse_jacobian[1] * r + cell_rule.inverse_jacobian[3] * s;
      }
    }

    const auto& local_faces = reference.faces;
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
      face_rule.minus = face_reference(group_index, face_index, quadrature_order, false);
      face_rule.weight_scale = 0.5 * cell_geometry.faces[face_index].measure;
      face_rule.lift_scale = face_rule.weight_scale / determinant;
      if (connection.is_boundary()) continue;

      const auto& permutation = connection.vertex_permutation;
      if (permutation.size() != 2 || permutation[0] > 1 || permutation[1] != 1 - permutation[0]) {
        throw std::logic_error("DGIntegration expects a two-vertex face permutation");
      }
      face_rule.plus = face_reference(neighbor_group, connection.neighbor_face, quadrature_order,
                                      permutation[0] == 1);
    }
  }
}

DGIntegration::Matrix DGIntegration::CellRule::mass() const {
  if (!reference) return {};
  return volume_jacobian * reference->mass;
}

DGIntegration::Matrix DGIntegration::CellRule::inverse_mass() const {
  if (!reference) return {};
  return reference->inverse_mass / volume_jacobian;
}

DGIntegration::Matrix DGIntegration::CellRule::stiffness_x() const {
  if (!reference) return {};
  return volume_jacobian * (inverse_jacobian[0] * reference->stiffness_r +
                            inverse_jacobian[2] * reference->stiffness_s);
}

DGIntegration::Matrix DGIntegration::CellRule::stiffness_y() const {
  if (!reference) return {};
  return volume_jacobian * (inverse_jacobian[1] * reference->stiffness_r +
                            inverse_jacobian[3] * reference->stiffness_s);
}

const DGIntegration::CellRule& DGIntegration::cell(std::size_t mesh_cell_index) const {
  return m_cells.at(mesh_cell_index);
}

}  // namespace oiseau::dg
