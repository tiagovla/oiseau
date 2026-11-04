// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/solver/linear_dg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "oiseau/dg/dg_connectivity.hpp"
#include "oiseau/dg/dg_geometry.hpp"
#include "oiseau/dg/dg_integration.hpp"

namespace oiseau::dg::solver {
namespace {

using FluxMatrices = std::array<std::vector<std::vector<double>>, 2>;

}  // namespace

struct LinearDG2D::Impl {
  struct FaceFluxOperators {
    std::vector<double> normal_flux_matrix;  // [equation, field]
  };

  struct CellFluxOperators {
    std::vector<FaceFluxOperators> faces;
  };

  std::shared_ptr<DGSpace> owned_space;
  const DGSpace* space = nullptr;
  DGConnectivity connectivity;
  DGGeometry geometry;
  DGIntegration integration;
  EquationSystem equations;
  NumericalFlux numerical_flux;
  BoundaryState boundary_state;
  FluxMatrices flux_matrices;
  std::vector<CellFluxOperators> normal_flux_operators;

  Impl(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem input_equations,
       NumericalFlux input_numerical_flux, BoundaryState input_boundary_state)
      : owned_space(std::make_shared<DGSpace>(mesh, orders)),
        space(owned_space.get()),
        connectivity(mesh),
        geometry(mesh),
        integration(mesh, *space, connectivity, geometry),
        equations(std::move(input_equations)),
        numerical_flux(std::move(input_numerical_flux)),
        boundary_state(std::move(input_boundary_state)) {
    initialize_flux_matrices();
    initialize_normal_flux_operators();
  }

  Impl(mesh::Mesh& mesh, const DGSpace& input_space, EquationSystem input_equations,
       NumericalFlux input_numerical_flux, BoundaryState input_boundary_state)
      : space(&input_space),
        connectivity(mesh),
        geometry(mesh),
        integration(mesh, input_space, connectivity, geometry),
        equations(std::move(input_equations)),
        numerical_flux(std::move(input_numerical_flux)),
        boundary_state(std::move(input_boundary_state)) {
    initialize_flux_matrices();
    initialize_normal_flux_operators();
  }

  void initialize_flux_matrices() {
    if (!numerical_flux) {
      throw std::invalid_argument("LinearDG2D requires a numerical flux policy");
    }

    const auto fields = equations.fields();
    if (std::ranges::any_of(equations.equations(), [](const Equation& equation) {
          return std::ranges::any_of(equation.terms, [](const LinearTerm& term) {
            return term.derivative == Derivative::Z && term.coefficient != 0.0;
          });
        })) {
      throw std::invalid_argument("LinearDG2D only supports x and y derivatives");
    }
    for (auto& matrix : flux_matrices) {
      matrix.assign(fields.size(), std::vector<double>(fields.size(), 0.0));
    }
    const std::array<Derivative, 2> directions = {Derivative::X, Derivative::Y};
    for (std::size_t direction = 0; direction < directions.size(); ++direction) {
      for (std::size_t equation = 0; equation < fields.size(); ++equation) {
        for (std::size_t field = 0; field < fields.size(); ++field) {
          flux_matrices[direction][equation][field] =
              equations.flux_coefficient(fields[equation], directions[direction], fields[field]);
        }
      }
    }
  }

  void initialize_normal_flux_operators() {
    const auto fields = equations.fields().size();
    normal_flux_operators.reserve(integration.cells().size());
    for (const auto& cell_rule : integration.cells()) {
      auto& cell = normal_flux_operators.emplace_back();
      if (cell_rule.group_index == DGConnectivity::no_neighbor) continue;

      cell.faces.reserve(cell_rule.faces.size());
      for (const auto& face_rule : cell_rule.faces) {
        auto& face = cell.faces.emplace_back();
        face.normal_flux_matrix.resize(fields * fields);
        for (std::size_t equation = 0; equation < fields; ++equation) {
          for (std::size_t field = 0; field < fields; ++field) {
            face.normal_flux_matrix[equation * fields + field] =
                face_rule.outward_normal[0] * flux_matrices[0][equation][field] +
                face_rule.outward_normal[1] * flux_matrices[1][equation][field];
          }
        }
      }
    }
  }

  void validate_state(const State& state) const {
    const auto groups = space->element_groups();
    if (state.groups.size() != groups.size()) {
      throw std::invalid_argument("DGState group count does not match DGSpace");
    }
    for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
      const auto& values = state.groups[group_index];
      const auto& nodes = groups[group_index].nodes;
      if (values.dimension() != 3 || values.shape()[0] != nodes.shape()[0] ||
          values.shape()[1] != nodes.shape()[1] || values.shape()[2] != equations.fields().size()) {
        throw std::invalid_argument(
            "DGState arrays must have shape [group cells, group nodes, number of fields]");
      }
      const auto& strides = values.strides();
      if ((values.shape()[2] > 1 && strides[2] != 1) ||
          (values.shape()[1] > 1 && std::cmp_not_equal(strides[1], values.shape()[2])) ||
          (values.shape()[0] > 1 && strides[0] != values.shape()[1] * values.shape()[2])) {
        throw std::invalid_argument("DGState arrays must use contiguous row-major storage");
      }
    }
  }

  State zero_state() const {
    State state;
    const auto fields = equations.fields().size();
    const auto groups = space->element_groups();
    state.groups.reserve(groups.size());
    for (const auto& group : groups) {
      const auto& shape = group.nodes.shape();
      const std::array<std::size_t, 3> state_shape = {shape[0], shape[1], fields};
      state.groups.emplace_back(xt::zeros<double>(state_shape));
    }
    return state;
  }

  void set_initial_condition(State& state, const Field& field,
                             const InitialCondition& value) const {
    if (!value) throw std::invalid_argument("Initial condition callback cannot be empty");
    validate_state(state);
    const auto field_index = equations.field_index(field);
    const auto groups = space->element_groups();
    for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
      const auto& nodes = groups[group_index].nodes;
      for (std::size_t cell = 0; cell < nodes.shape()[0]; ++cell) {
        for (std::size_t node = 0; node < nodes.shape()[1]; ++node) {
          state.groups[group_index](cell, node, field_index) =
              value(nodes(cell, node, 0), nodes(cell, node, 1));
        }
      }
    }
  }

  State initial_state(const Field& field, const InitialCondition& value) const {
    auto state = zero_state();
    set_initial_condition(state, field, value);
    return state;
  }

  State rhs(const State& state) const {
    validate_state(state);
    auto output = zero_state();
    const auto fields = equations.fields().size();
    std::vector<double> left(fields), right(fields), left_flux(fields), right_flux(fields);
    std::vector<double> volume_flux_x, volume_flux_y;

    const auto cell_rules = integration.cells();
    for (std::size_t cell_index = 0; cell_index < cell_rules.size(); ++cell_index) {
      const auto& cell_rule = cell_rules[cell_index];
      if (cell_rule.group_index == DGConnectivity::no_neighbor) continue;
      const auto& operators = normal_flux_operators[cell_index];
      const auto group = cell_rule.group_index;
      const auto row = cell_rule.group_cell_index;
      const auto& cell_state = state.groups[group];
      const auto node_count = cell_rule.node_count();
      const auto* cell_values = cell_state.data() + row * node_count * fields;
      auto* cell_derivative = output.groups[group].data() + row * node_count * fields;

      volume_flux_x.resize(node_count * fields);
      volume_flux_y.resize(node_count * fields);
      for (std::size_t trial = 0; trial < node_count; ++trial) {
        for (std::size_t equation = 0; equation < fields; ++equation) {
          double x = 0.0;
          double y = 0.0;
          for (std::size_t field = 0; field < fields; ++field) {
            const auto value = cell_values[trial * fields + field];
            x += flux_matrices[0][equation][field] * value;
            y += flux_matrices[1][equation][field] * value;
          }
          volume_flux_x[trial * fields + equation] = x;
          volume_flux_y[trial * fields + equation] = y;
        }
      }
      for (std::size_t test = 0; test < node_count; ++test) {
        for (std::size_t equation = 0; equation < fields; ++equation) {
          double value = 0.0;
          for (std::size_t trial = 0; trial < node_count; ++trial) {
            const auto index = test * node_count + trial;
            value += cell_rule.weak_derivative_x[index] * volume_flux_x[trial * fields + equation] +
                     cell_rule.weak_derivative_y[index] * volume_flux_y[trial * fields + equation];
          }
          cell_derivative[test * fields + equation] = value;
        }
      }

      for (std::size_t face_index = 0; face_index < cell_rule.faces.size(); ++face_index) {
        const auto& face = cell_rule.faces[face_index];
        const auto& face_operators = operators.faces[face_index];
        const auto* neighbor_rule =
            face.is_boundary() ? &cell_rule : &integration.cell(face.neighbor_cell);
        const auto neighbor_group = neighbor_rule->group_index;
        const auto neighbor_row = neighbor_rule->group_cell_index;
        const auto& neighbor_state = state.groups[neighbor_group];
        const auto neighbor_nodes = neighbor_state.shape()[1];
        const auto* neighbor_values =
            neighbor_state.data() + neighbor_row * neighbor_nodes * fields;
        const auto face_points = face.point_count();

        for (std::size_t q = 0; q < face_points; ++q) {
          std::ranges::fill(left, 0.0);
          std::ranges::fill(right, 0.0);
          const auto* minus_trace = face.minus->trace.data() + q * node_count;
          for (std::size_t node = 0; node < node_count; ++node) {
            const auto basis = minus_trace[node];
            for (std::size_t field = 0; field < fields; ++field) {
              left[field] += basis * cell_values[node * fields + field];
            }
          }

          if (face.is_boundary()) {
            right = boundary_state ? boundary_state(left, face.outward_normal) : left;
            if (right.size() != fields) {
              throw std::invalid_argument(
                  "Boundary state must return one value for every equation field");
            }
          } else {
            const auto* plus_trace = face.plus->trace.data() + q * neighbor_nodes;
            for (std::size_t node = 0; node < neighbor_nodes; ++node) {
              const auto basis = plus_trace[node];
              for (std::size_t field = 0; field < fields; ++field) {
                right[field] += basis * neighbor_values[node * fields + field];
              }
            }
          }

          for (std::size_t equation = 0; equation < fields; ++equation) {
            const auto* normal_matrix =
                face_operators.normal_flux_matrix.data() + equation * fields;
            double interior_flux = 0.0;
            double exterior_flux = 0.0;
            for (std::size_t field = 0; field < fields; ++field) {
              interior_flux += normal_matrix[field] * left[field];
              exterior_flux += normal_matrix[field] * right[field];
            }
            left_flux[equation] = interior_flux;
            right_flux[equation] = exterior_flux;
          }
          const auto flux = numerical_flux(left, right, face.outward_normal, left_flux, right_flux);
          if (flux.size() != fields) {
            throw std::invalid_argument("Numerical flux must return one value per equation field");
          }
          const auto* lift = face.minus->lift.data() + q * node_count;
          for (std::size_t test = 0; test < node_count; ++test) {
            for (std::size_t equation = 0; equation < fields; ++equation) {
              cell_derivative[test * fields + equation] -=
                  face.lift_scale * lift[test] * flux[equation];
            }
          }
        }
      }
    }
    return output;
  }
};

LinearDG2D::LinearDG2D(mesh::Mesh& mesh, const std::vector<unsigned>& orders,
                       EquationSystem equations, WaveSpeed wave_speed, BoundaryState boundary_state)
    : LinearDG2D(mesh, orders, std::move(equations), rusanov(std::move(wave_speed)),
                 std::move(boundary_state)) {}

LinearDG2D::LinearDG2D(mesh::Mesh& mesh, const std::vector<unsigned>& orders,
                       EquationSystem equations, NumericalFlux numerical_flux,
                       BoundaryState boundary_state)
    : m_impl(std::make_unique<Impl>(mesh, orders, std::move(equations), std::move(numerical_flux),
                                    std::move(boundary_state))) {}

LinearDG2D::LinearDG2D(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
                       WaveSpeed wave_speed, BoundaryState boundary_state)
    : LinearDG2D(mesh, space, std::move(equations), rusanov(std::move(wave_speed)),
                 std::move(boundary_state)) {}

LinearDG2D::LinearDG2D(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
                       NumericalFlux numerical_flux, BoundaryState boundary_state)
    : m_impl(std::make_unique<Impl>(mesh, space, std::move(equations), std::move(numerical_flux),
                                    std::move(boundary_state))) {}

LinearDG2D::~LinearDG2D() = default;
LinearDG2D::LinearDG2D(LinearDG2D&&) noexcept = default;
LinearDG2D& LinearDG2D::operator=(LinearDG2D&&) noexcept = default;

const DGSpace& LinearDG2D::space() const { return *m_impl->space; }
const EquationSystem& LinearDG2D::equations() const { return m_impl->equations; }
LinearDG2D::NumericalFlux LinearDG2D::rusanov(WaveSpeed wave_speed) {
  if (!wave_speed) throw std::invalid_argument("Rusanov flux requires a wave-speed bound");
  return [wave_speed = std::move(wave_speed)](
             std::span<const double> left, std::span<const double> right, const Normal& normal,
             std::span<const double> left_flux, std::span<const double> right_flux) {
    if (left.size() != right.size() || left.size() != left_flux.size() ||
        left.size() != right_flux.size()) {
      throw std::invalid_argument("Rusanov flux state and physical flux sizes must match");
    }
    const auto speed = wave_speed(normal);
    if (!std::isfinite(speed) || speed < 0.0) {
      throw std::invalid_argument("Rusanov wave-speed bound must be finite and nonnegative");
    }
    std::vector<double> flux(left.size());
    for (std::size_t field = 0; field < flux.size(); ++field) {
      flux[field] =
          0.5 * (left_flux[field] + right_flux[field]) - 0.5 * speed * (right[field] - left[field]);
    }
    return flux;
  };
}
LinearDG2D::State LinearDG2D::zero_state() const { return m_impl->zero_state(); }
LinearDG2D::State LinearDG2D::initial_state(const Field& field,
                                            const InitialCondition& value) const {
  return m_impl->initial_state(field, value);
}
void LinearDG2D::set_initial_condition(State& state, const Field& field,
                                       const InitialCondition& value) const {
  m_impl->set_initial_condition(state, field, value);
}
LinearDG2D::State LinearDG2D::rhs(const State& state) const { return m_impl->rhs(state); }

}  // namespace oiseau::dg::solver
