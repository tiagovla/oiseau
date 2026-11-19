// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_integration.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
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
#include <xtensor/misc/xmanipulation.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/dg/detail/cell_mapping.hpp"
#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/dg_space.hpp"
#include "oiseau/dg/nodal/ref_element.hpp"
#include "oiseau/dg/ref_element_type.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/utils/integration.hpp"

namespace oiseau::dg {
namespace detail {

using Matrix = DGIntegration::Matrix;
using Point = std::array<double, 3>;
using FacePermutation = std::array<std::size_t, 4>;

constexpr FacePermutation identity_permutation{0, 1, 2, 3};
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
  const auto& vandermonde = element->v();
  auto inverse_vandermonde = xt::linalg::inv(vandermonde);

  auto volume = std::make_shared<DGIntegration::ReferenceRule>();
  volume->mass = xt::linalg::dot(xt::transpose(inverse_vandermonde), inverse_vandermonde);
  volume->inverse_mass = xt::linalg::dot(vandermonde, xt::transpose(vandermonde));
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
          .inverse_vandermonde = std::move(inverse_vandermonde)};
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

DGIntegration::QuadratureRule make_face_quadrature(std::size_t corners, unsigned order) {
  if (corners == 1) {
    return {.nodes = Matrix::from_shape({1, 0}), .weights = Matrix{1.0}};
  }
  if (order > (static_cast<unsigned>(std::numeric_limits<int>::max()) - 2U) / 2U) {
    throw std::invalid_argument("Polynomial order is too large for DG face quadrature");
  }
  const auto degree = static_cast<int>(2U * order + 1U);
  if (corners == 3) {
    auto [nodes, weights] = utils::integration::cubature(degree);
    return {.nodes = std::move(nodes), .weights = std::move(weights)};
  }

  auto [nodes, weights] = utils::integration::quadrature(degree);
  if (corners == 2) return {.nodes = std::move(nodes), .weights = std::move(weights)};
  if (corners != 4) throw std::logic_error("DGIntegration does not support this face shape");

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

Matrix face_points(const ReferenceOperators& reference, const std::vector<int>& vertices,
                   const FacePermutation& permutation,
                   const DGIntegration::QuadratureRule& quadrature) {
  if (vertices.size() == 1) return Matrix{reference.vertices[vertices.front()][0]};
  const auto dimension = reference.element->d().shape()[2];
  Matrix points = xt::zeros<double>({quadrature.weights.size(), dimension});
  for (std::size_t q = 0; q < quadrature.weights.size(); ++q) {
    const auto weights = corner_weights(vertices.size(), quadrature.nodes, q);
    for (std::size_t corner = 0; corner < vertices.size(); ++corner) {
      const auto& vertex = reference.vertices[vertices[permutation[corner]]];
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
  const auto& inverse_jacobian = cell.spatial_inverse_jacobian;
  Matrix stiffness = inverse_jacobian[d] * reference.stiffness_r;
  if (has_y) stiffness += inverse_jacobian[3 + d] * reference.stiffness_s;
  if (has_z) stiffness += inverse_jacobian[6 + d] * reference.stiffness_t;
  stiffness *= cell.volume_jacobian;
  return stiffness;
}

}  // namespace detail

using namespace detail;

DGIntegration::DGIntegration(const mesh::Mesh& mesh, const DGSpace& space,
                             const DGConnectivity& connectivity, const DGGeometry& geometry) {
  if (&space.mesh() != &mesh || &connectivity.mesh() != &mesh || &geometry.mesh() != &mesh) {
    throw std::invalid_argument("DGIntegration inputs must refer to the same mesh");
  }

  const auto groups = space.element_groups();
  const auto cell_types = mesh.topology().cell_types();
  const auto connectivity_data = mesh.topology().conn();
  const auto cell_count = mesh.topology().n_cells();
  const auto face_dimension = mesh.topology().dimension() - 1;

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
    reference.faces = cell_types[first_cell]->get_entity_vertices(face_dimension);
    for (std::size_t row = 0; row < group.cell_indices.size(); ++row) {
      const auto cell_index = group.cell_indices[row];
      cell_group[cell_index] = group_index;
      cell_group_row[cell_index] = row;
    }
  }

  std::map<std::pair<std::size_t, unsigned>, std::shared_ptr<const QuadratureRule>> quadratures;
  using FaceKey = std::tuple<std::size_t, std::size_t, unsigned, FacePermutation>;
  std::map<FaceKey, std::shared_ptr<const ReferenceFaceRule>> face_rules;
  const auto face_reference = [&](std::size_t group_index, std::size_t face_index, unsigned order,
                                  const FacePermutation& permutation) {
    const auto& reference = reference_operators[group_index];
    const auto& vertices = reference.faces.at(face_index);
    const auto face_order = vertices.size() == 1 ? 0U : order;
    const FaceKey key{group_index, face_index, face_order, permutation};
    if (const auto found = face_rules.find(key); found != face_rules.end()) return found->second;

    auto& quadrature = quadratures[{vertices.size(), face_order}];
    if (!quadrature) {
      quadrature =
          std::make_shared<const QuadratureRule>(make_face_quadrature(vertices.size(), face_order));
    }
    const auto points = face_points(reference, vertices, permutation, *quadrature);
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
    cell_rule.spatial_inverse_jacobian = cell_geometry.inverse_jacobian;

    const auto nodes = cell_rule.node_count();
    const auto& inverse_jacobian = cell_rule.spatial_inverse_jacobian;
    const auto has_y = cell_geometry.dimension >= 2;
    const auto has_z = cell_geometry.dimension == 3;
    const auto set_weak_derivative = [&](std::size_t d, std::vector<double>& derivative) {
      derivative.resize(nodes * nodes);
      for (std::size_t test = 0; test < nodes; ++test) {
        for (std::size_t trial = 0; trial < nodes; ++trial) {
          auto entry = inverse_jacobian[d] * reference.volume->weak_derivative_r(test, trial);
          if (has_y) {
            entry += inverse_jacobian[3 + d] * reference.volume->weak_derivative_s(test, trial);
          }
          if (has_z) {
            entry += inverse_jacobian[6 + d] * reference.volume->weak_derivative_t(test, trial);
          }
          derivative[test * nodes + trial] = entry;
        }
      }
    };
    set_weak_derivative(0, cell_rule.weak_derivative_x);
    if (has_y) set_weak_derivative(1, cell_rule.weak_derivative_y);
    if (has_z) {
      set_weak_derivative(2, cell_rule.weak_derivative_z);
    }

    const auto& local_faces = reference.faces;
    const auto face_connections = connectivity.faces(cell_index);
    cell_rule.faces.reserve(local_faces.size());
    for (std::size_t face_index = 0; face_index < local_faces.size(); ++face_index) {
      const auto& connection = face_connections[face_index];
      auto& face_rule = cell_rule.faces.emplace_back();
      const auto& face_geometry = cell_geometry.faces[face_index];
      face_rule.outward_normal = {face_geometry.outward_normal[0], face_geometry.outward_normal[1]};
      face_rule.normal_z = face_geometry.outward_normal[2];
      face_rule.neighbor_cell = connection.neighbor_cell;

      const auto corners = local_faces[face_index].size();
      auto neighbor_group = group_index;
      auto quadrature_order = groups[group_index].order;
      auto minus_permutation = identity_permutation;
      auto plus_permutation = identity_permutation;
      if (!connection.is_boundary()) {
        neighbor_group = cell_group.at(connection.neighbor_cell);
        if (neighbor_group == DGConnectivity::no_neighbor) {
          throw std::invalid_argument("DGIntegration neighbor cell is absent from DGSpace");
        }
        quadrature_order = std::max(quadrature_order, groups[neighbor_group].order);
        const auto permutation = face_permutation(connection.vertex_permutation, corners);
        // Both directions use one 3D face orientation, even for nonsymmetric cubature.
        if (has_z && connection.neighbor_cell < cell_index) {
          for (std::size_t corner = 0; corner < corners; ++corner) {
            minus_permutation[permutation[corner]] = corner;
          }
        } else {
          plus_permutation = permutation;
        }
      }
      face_rule.minus =
          face_reference(group_index, face_index, quadrature_order, minus_permutation);
      face_rule.weight_scale = face_geometry.measure / canonical_face_measures[corners - 1];
      face_rule.lift_scale = face_rule.weight_scale / determinant;
      if (!connection.is_boundary()) {
        face_rule.plus = face_reference(neighbor_group, connection.neighbor_face, quadrature_order,
                                        plus_permutation);
      }
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
