// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/integration.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/connectivity.hpp"
#include "oiseau/dg/detail/cell_mapping.hpp"
#include "oiseau/dg/detail/integration.hpp"
#include "oiseau/dg/geometry.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/nodal/ref_element_type.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/utils/integration.hpp"
#include "oiseau/utils/logging.hpp"

namespace oiseau::dg {
namespace detail {

using Matrix = DGIntegration::Matrix;
using Point = std::array<double, 3>;

constexpr FacePermutation identity_permutation{0, 1, 2, 3};
constexpr std::array canonical_face_kinds{mesh::CellKind::Point, mesh::CellKind::Interval,
                                          mesh::CellKind::Triangle, mesh::CellKind::Quadrilateral};
constexpr std::array canonical_face_measures{1.0, 2.0, 2.0, 4.0};

struct ReferenceOperators {
  std::shared_ptr<nodal::RefElement> element;
  std::shared_ptr<const DGIntegration::ReferenceRule> volume;
  Matrix inverse_vandermonde;
  std::vector<Point> vertices;
  std::vector<std::vector<int>> faces;
};

ReferenceOperators make_reference_operators(mesh::CellKind kind, unsigned order) {
  auto element = nodal::get_ref_element(ref_type_for(kind), order);
  auto inverse_vandermonde = element->inv_v();

  auto volume = std::make_shared<DGIntegration::ReferenceRule>();
  volume->mass = element->mass();
  volume->inverse_mass = element->inv_mass();
  const auto& derivatives = element->d();
  const auto add_operators = [&](const auto& values, Matrix& derivative, Matrix& stiffness,
                                 Matrix& weak_derivative) {
    derivative = values;
    stiffness = xt::linalg::dot(xt::transpose(derivative), volume->mass);
    weak_derivative = xt::linalg::dot(volume->inverse_mass, stiffness);
  };
  if (derivatives.dimension() == 2) {
    add_operators(derivatives, volume->derivative_r, volume->stiffness_r,
                  volume->weak_derivative_r);
  } else {
    add_operators(xt::view(derivatives, xt::all(), xt::all(), 0), volume->derivative_r,
                  volume->stiffness_r, volume->weak_derivative_r);
    add_operators(xt::view(derivatives, xt::all(), xt::all(), 1), volume->derivative_s,
                  volume->stiffness_s, volume->weak_derivative_s);
    if (derivatives.shape()[2] == 3) {
      add_operators(xt::view(derivatives, xt::all(), xt::all(), 2), volume->derivative_t,
                    volume->stiffness_t, volume->weak_derivative_t);
    }
  }
  return {.element = std::move(element),
          .volume = std::move(volume),
          .inverse_vandermonde = std::move(inverse_vandermonde),
          .vertices = {},
          .faces = {}};
}

std::vector<Point> reference_vertex_coordinates(mesh::CellType type,
                                                const std::vector<std::size_t>& mesh_vertices) {
  const auto vertices = detail::reference_ordered_vertices(type, mesh_vertices);
  const auto reference = nodal::get_ref_element(ref_type_for(type->kind()), 1);
  const auto& nodes = reference->r();
  const auto dimension = static_cast<std::size_t>(type->dimension());
  std::vector<Point> output(mesh_vertices.size());
  for (std::size_t reference_vertex = 0; reference_vertex < vertices.size(); ++reference_vertex) {
    const auto position = std::ranges::find(mesh_vertices, vertices[reference_vertex]);
    if (position == mesh_vertices.end()) {
      throw std::logic_error("Reference vertex is missing from mesh connectivity");
    }
    const auto mesh_vertex = static_cast<std::size_t>(position - mesh_vertices.begin());
    for (std::size_t d = 0; d < dimension; ++d) {
      output[mesh_vertex][d] =
          dimension == 1 ? nodes(reference_vertex) : nodes(reference_vertex, d);
    }
  }
  return output;
}

DGIntegration::QuadratureRule make_face_quadrature(mesh::CellKind face_kind, unsigned order) {
  using mesh::CellKind;
  if (face_kind == CellKind::Point) {
    return {.nodes = Matrix::from_shape({1, 0}), .weights = Matrix{1.0}};
  }
  if (order > (static_cast<unsigned>(std::numeric_limits<int>::max()) - 2U) / 2U) {
    throw std::invalid_argument("Polynomial order is too large for DG face quadrature");
  }
  const auto degree = static_cast<int>(2U * order + 1U);
  if (face_kind == CellKind::Triangle) {
    auto [nodes, weights] = utils::integration::triangle_cubature(degree);
    return {.nodes = std::move(nodes), .weights = std::move(weights)};
  }

  auto [nodes, weights] = utils::integration::quadrature(degree);
  if (face_kind == CellKind::Interval) {
    return {.nodes = std::move(nodes), .weights = std::move(weights)};
  }
  if (face_kind != CellKind::Quadrilateral) {
    throw std::logic_error("DGIntegration does not support this face shape");
  }

  const auto count = nodes.size();
  Matrix square_nodes = Matrix::from_shape({count * count, std::size_t{2}});
  Matrix square_weights = Matrix::from_shape({count * count});
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = 0; j < count; ++j) {
      const auto q = i * count + j;
      square_nodes(q, 0) = nodes(i);
      square_nodes(q, 1) = nodes(j);
      square_weights(q) = weights(i) * weights(j);
    }
  }
  return {.nodes = std::move(square_nodes), .weights = std::move(square_weights)};
}

std::array<double, 4> corner_weights(std::size_t corners, const Matrix& nodes, std::size_t q) {
  if (corners == 2) {
    const auto u = 0.5 * (1.0 + nodes(q));
    return {1.0 - u, u, 0.0, 0.0};
  }
  const auto u = 0.5 * (1.0 + nodes(q, 0));
  const auto v = 0.5 * (1.0 + nodes(q, 1));
  if (corners == 3) return {1.0 - u - v, u, v, 0.0};
  return {(1.0 - u) * (1.0 - v), u * (1.0 - v), u * v, (1.0 - u) * v};
}

Matrix face_points(std::span<const Point> cell_vertices, std::span<const int> face_vertices,
                   std::size_t dimension, const FacePermutation& permutation,
                   const DGIntegration::QuadratureRule& quadrature) {
  if (face_vertices.size() == 1) return Matrix{cell_vertices[face_vertices.front()][0]};
  Matrix points = xt::zeros<double>({quadrature.weights.size(), dimension});
  for (std::size_t q = 0; q < quadrature.weights.size(); ++q) {
    const auto weights = corner_weights(face_vertices.size(), quadrature.nodes, q);
    for (std::size_t corner = 0; corner < face_vertices.size(); ++corner) {
      const auto& vertex = cell_vertices[face_vertices[permutation[corner]]];
      for (std::size_t d = 0; d < dimension; ++d) {
        points(q, d) += weights[corner] * vertex[d];
      }
    }
  }
  return points;
}

FacePermutation face_permutation(std::span<const std::size_t> neighbor_positions,
                                 std::size_t corners) {
  if (neighbor_positions.size() != corners) {
    throw std::logic_error("DGIntegration face permutation has the wrong size");
  }

  FacePermutation permutation = identity_permutation;
  std::array<bool, 4> seen{};
  for (std::size_t i = 0; i < corners; ++i) {
    const auto position = neighbor_positions[i];
    if (position >= corners || seen[position]) {
      throw std::logic_error("DGIntegration face permutation is invalid");
    }
    permutation[i] = position;
    seen[position] = true;
  }
  return permutation;
}

Matrix interpolation_on_points(const ReferenceOperators& operators, const Matrix& points) {
  const auto vandermonde = operators.element->vandermonde(points);
  return xt::linalg::dot(vandermonde, operators.inverse_vandermonde);
}

Matrix physical_stiffness(const DGIntegration::CellRule& cell, std::size_t d) {
  if (!cell.reference) return {};
  const auto& reference = *cell.reference;
  const auto has_y = reference.stiffness_s.dimension() == 2;
  const auto has_z = reference.stiffness_t.dimension() == 2;
  if ((d == 1 && !has_y) || (d == 2 && !has_z)) return {};
  const auto& inverse_jacobian = cell.inverse_jacobian;
  Matrix stiffness = inverse_jacobian[d] * reference.stiffness_r;
  if (has_y) stiffness += inverse_jacobian[3 + d] * reference.stiffness_s;
  if (has_z) stiffness += inverse_jacobian[6 + d] * reference.stiffness_t;
  stiffness *= cell.volume_jacobian;
  return stiffness;
}

}  // namespace detail

using namespace detail;

namespace {

std::vector<ReferenceOperators> make_group_operators(const mesh::Topology& topology,
                                                     std::span<const DGSpace::ElementGroup> groups,
                                                     int face_dimension,
                                                     std::span<DGIntegration::CellRule> cells) {
  const auto cell_types = topology.cell_types();
  const auto connectivity = topology.conn();
  std::vector<ReferenceOperators> reference_operators;
  reference_operators.reserve(groups.size());
  for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
    const auto& group = groups[group_index];
    auto& reference =
        reference_operators.emplace_back(make_reference_operators(group.kind, group.order));
    const auto first_cell = group.cell_indices.front();
    reference.vertices =
        reference_vertex_coordinates(cell_types[first_cell], connectivity[first_cell]);
    reference.faces = cell_types[first_cell]->get_entity_vertices(face_dimension);
    for (std::size_t row = 0; row < group.cell_indices.size(); ++row) {
      auto& cell_rule = cells[group.cell_indices[row]];
      cell_rule.group_index = group_index;
      cell_rule.group_cell_index = row;
    }
  }
  return reference_operators;
}

void set_weak_derivative(std::vector<double>& derivative, const DGIntegration::CellRule& cell,
                         std::size_t dimension, std::size_t d) {
  const auto nodes = cell.node_count();
  const auto& reference = *cell.reference;
  const auto& inverse_jacobian = cell.inverse_jacobian;
  const auto has_y = dimension >= 2;
  const auto has_z = dimension == 3;
  derivative.resize(nodes * nodes);
  for (std::size_t test = 0; test < nodes; ++test) {
    for (std::size_t trial = 0; trial < nodes; ++trial) {
      auto entry = inverse_jacobian[d] * reference.weak_derivative_r(test, trial);
      if (has_y) {
        entry += inverse_jacobian[3 + d] * reference.weak_derivative_s(test, trial);
      }
      if (has_z) {
        entry += inverse_jacobian[6 + d] * reference.weak_derivative_t(test, trial);
      }
      derivative[test * nodes + trial] = entry;
    }
  }
}

void set_volume_rule(DGIntegration::CellRule& cell,
                     const std::shared_ptr<const DGIntegration::ReferenceRule>& reference,
                     const DGGeometry::CellGeometry& geometry) {
  cell.reference = reference;
  cell.volume_jacobian = std::abs(geometry.determinant);
  cell.inverse_jacobian = geometry.inverse_jacobian;

  const auto dimension = geometry.dimension;
  set_weak_derivative(cell.weak_derivative_x, cell, dimension, 0);
  if (dimension >= 2) set_weak_derivative(cell.weak_derivative_y, cell, dimension, 1);
  if (dimension == 3) set_weak_derivative(cell.weak_derivative_z, cell, dimension, 2);
}

FacePermutation inverse_face_permutation(const FacePermutation& permutation, std::size_t corners) {
  auto inverse = identity_permutation;
  for (std::size_t corner = 0; corner < corners; ++corner) {
    inverse[permutation[corner]] = corner;
  }
  return inverse;
}

class FaceRuleBuilder {
  using QuadratureRule = DGIntegration::QuadratureRule;
  using ReferenceFaceRule = DGIntegration::ReferenceFaceRule;
  using FaceKey = std::tuple<std::size_t, std::size_t, unsigned, FacePermutation>;

 public:
  FaceRuleBuilder(std::span<const DGSpace::ElementGroup> groups,
                  std::span<const ReferenceOperators> references, std::size_t dimension)
      : m_groups(groups), m_references(references), m_dimension(dimension) {}

  void add_cell_faces(std::vector<DGIntegration::CellRule>& cells, std::size_t cell_index,
                      std::span<const DGConnectivity::FaceConnection> face_connections,
                      const DGGeometry::CellGeometry& geometry) {
    auto& cell_rule = cells[cell_index];
    const auto group_index = cell_rule.group_index;
    const auto& local_faces = m_references[group_index].faces;
    cell_rule.faces.reserve(local_faces.size());
    for (std::size_t face_index = 0; face_index < local_faces.size(); ++face_index) {
      const auto& connection = face_connections[face_index];
      auto& face_rule = cell_rule.faces.emplace_back();
      const auto& face_geometry = geometry.faces[face_index];
      face_rule.outward_normal = face_geometry.outward_normal;
      face_rule.neighbor_cell = connection.neighbor_cell;

      const auto corners = local_faces[face_index].size();
      auto neighbor_group = group_index;
      auto quadrature_order = m_groups[group_index].order;
      auto minus_permutation = identity_permutation;
      auto plus_permutation = identity_permutation;
      if (!connection.is_boundary()) {
        const auto& neighbor = cells.at(connection.neighbor_cell);
        if (!neighbor.has_group()) {
          throw std::invalid_argument("DGIntegration neighbor cell is absent from DGSpace");
        }
        neighbor_group = neighbor.group_index;
        quadrature_order = std::max(quadrature_order, m_groups[neighbor_group].order);
        const auto permutation = face_permutation(connection.vertex_permutation, corners);
        // Both directions use one 3D face orientation, even for nonsymmetric cubature.
        if (geometry.dimension == 3 && connection.neighbor_cell < cell_index) {
          minus_permutation = inverse_face_permutation(permutation, corners);
        } else {
          plus_permutation = permutation;
        }
      }
      face_rule.minus =
          face_reference(group_index, face_index, quadrature_order, minus_permutation);
      face_rule.weight_scale = face_geometry.measure / canonical_face_measures[corners - 1];
      face_rule.lift_scale = face_rule.weight_scale / cell_rule.volume_jacobian;
      if (!connection.is_boundary()) {
        face_rule.plus = face_reference(neighbor_group, connection.neighbor_face, quadrature_order,
                                        plus_permutation);
      }
    }
  }

 private:
  std::shared_ptr<const ReferenceFaceRule> face_reference(std::size_t group_index,
                                                          std::size_t face_index, unsigned order,
                                                          const FacePermutation& permutation) {
    const auto& reference = m_references[group_index];
    const auto& vertices = reference.faces.at(face_index);
    const auto face_kind = canonical_face_kinds.at(vertices.size() - 1);
    const auto face_order = face_kind == mesh::CellKind::Point ? 0U : order;
    const FaceKey key{group_index, face_index, face_order, permutation};
    if (const auto found = m_face_rules.find(key); found != m_face_rules.end()) {
      return found->second;
    }

    auto& quadrature = m_quadratures[{face_kind, face_order}];
    if (!quadrature) {
      quadrature =
          std::make_shared<const QuadratureRule>(make_face_quadrature(face_kind, face_order));
    }
    const auto points =
        face_points(reference.vertices, vertices, m_dimension, permutation, *quadrature);
    Matrix trace = interpolation_on_points(reference, points);
    Matrix lift = xt::linalg::dot(trace, xt::transpose(reference.volume->inverse_mass));
    for (std::size_t q = 0; q < quadrature->weights.size(); ++q) {
      for (std::size_t node = 0; node < lift.shape()[1]; ++node) {
        lift(q, node) *= quadrature->weights(q);
      }
    }
    auto rule = std::make_shared<const ReferenceFaceRule>(ReferenceFaceRule{
        .quadrature = quadrature, .trace = std::move(trace), .lift = std::move(lift)});
    m_face_rules.emplace(key, rule);
    return rule;
  }

  std::span<const DGSpace::ElementGroup> m_groups;
  std::span<const ReferenceOperators> m_references;
  std::size_t m_dimension;
  std::map<std::pair<mesh::CellKind, unsigned>, std::shared_ptr<const QuadratureRule>>
      m_quadratures;
  std::map<FaceKey, std::shared_ptr<const ReferenceFaceRule>> m_face_rules;
};

}  // namespace

DGIntegration::DGIntegration(const mesh::Mesh& mesh, const DGSpace& space,
                             const DGConnectivity& connectivity, const DGGeometry& geometry) {
  if (&space.mesh() != &mesh || &connectivity.mesh() != &mesh || &geometry.mesh() != &mesh) {
    throw std::invalid_argument("DGIntegration inputs must refer to the same mesh");
  }

  const auto groups = space.element_groups();
  const auto& topology = mesh.topology();
  const auto cell_count = topology.n_cells();
  const auto dimension = topology.dimension();
  LOG_DEBUG("Precomputing DG volume and face integration rules for {} groups", groups.size());
  m_cells.resize(cell_count);
  const auto reference_operators = make_group_operators(topology, groups, dimension - 1, m_cells);

  FaceRuleBuilder face_rules(groups, reference_operators, static_cast<std::size_t>(dimension));
  for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
    auto& cell_rule = m_cells[cell_index];
    if (!cell_rule.has_group()) continue;

    const auto& cell_geometry = geometry.cell(cell_index);
    set_volume_rule(cell_rule, reference_operators[cell_rule.group_index].volume, cell_geometry);
    face_rules.add_cell_faces(m_cells, cell_index, connectivity.faces(cell_index), cell_geometry);
  }
}

DGIntegration::DGIntegration(mesh::Mesh& mesh, const DGSpace& space)
    : DGIntegration(mesh, space, DGConnectivity(mesh), DGGeometry(mesh)) {}

DGIntegration::Matrix DGIntegration::CellRule::mass() const {
  if (!reference) return {};
  return volume_jacobian * reference->mass;
}

DGIntegration::Matrix DGIntegration::CellRule::inverse_mass() const {
  if (!reference) return {};
  return reference->inverse_mass / volume_jacobian;
}

DGIntegration::Matrix DGIntegration::CellRule::stiffness_x() const {
  return physical_stiffness(*this, 0);
}

DGIntegration::Matrix DGIntegration::CellRule::stiffness_y() const {
  return physical_stiffness(*this, 1);
}

DGIntegration::Matrix DGIntegration::CellRule::stiffness_z() const {
  return physical_stiffness(*this, 2);
}

const DGIntegration::CellRule& DGIntegration::cell(std::size_t mesh_cell_index) const {
  return m_cells.at(mesh_cell_index);
}

}  // namespace oiseau::dg
