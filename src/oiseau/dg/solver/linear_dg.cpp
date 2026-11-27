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
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "oiseau/dg/connectivity.hpp"
#include "oiseau/dg/equation.hpp"
#include "oiseau/dg/field.hpp"
#include "oiseau/dg/geometry.hpp"
#include "oiseau/dg/integration.hpp"
#include "oiseau/dg/space.hpp"
#include "oiseau/mesh/mesh.hpp"
#include "oiseau/mesh/topology.hpp"
#include "oiseau/utils/logging.hpp"

namespace oiseau::dg::solver {
template <std::size_t Dimension>
struct LinearDG<Dimension>::Impl {
  struct FaceFluxOperators {
    Normal normal{};
    std::vector<double> normal_flux_matrix;  // [equation, field]
  };

  struct CellFluxOperators {
    std::vector<FaceFluxOperators> faces;
  };

  struct NodalValues {
    const double* data;
    std::size_t node_count;
  };

  struct FaceScratch {
    explicit FaceScratch(std::size_t fields)
        : left(fields), right(fields), left_flux(fields), right_flux(fields) {}

    std::vector<double> left;
    std::vector<double> right;
    std::vector<double> left_flux;
    std::vector<double> right_flux;
  };

  using ReferenceFlux = std::array<std::vector<double>, Dimension>;

  std::shared_ptr<DGSpace> owned_space;
  const DGSpace* space = nullptr;
  DGConnectivity connectivity;
  DGGeometry geometry;
  DGIntegration integration;
  EquationSystem equations;
  NumericalFlux numerical_flux;
  BoundaryState boundary_state;
  std::array<std::vector<std::vector<double>>, Dimension> flux_matrices;
  std::vector<CellFluxOperators> normal_flux_operators;

  Impl(mesh::Mesh& mesh, const std::vector<unsigned>& orders, EquationSystem input_equations,
       NumericalFlux input_numerical_flux, BoundaryState input_boundary_state)
      : owned_space(std::make_shared<DGSpace>(check_mesh_dimension(mesh), orders)),
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
        connectivity(check_mesh_dimension(mesh)),
        geometry(mesh),
        integration(mesh, input_space, connectivity, geometry),
        equations(std::move(input_equations)),
        numerical_flux(std::move(input_numerical_flux)),
        boundary_state(std::move(input_boundary_state)) {
    initialize_flux_matrices();
    initialize_normal_flux_operators();
  }

  static mesh::Mesh& check_mesh_dimension(mesh::Mesh& mesh) {
    if (mesh.topology().dimension() != static_cast<int>(Dimension)) {
      throw std::invalid_argument("LinearDG mesh dimension must match the solver dimension");
    }
    return mesh;
  }

  void validate_flux_policy() const {
    if (!numerical_flux) {
      throw std::invalid_argument("LinearDG requires a numerical flux policy");
    }

    if constexpr (Dimension < 3) {
      if (std::ranges::any_of(equations.equations(), [](const Equation& equation) {
            return std::ranges::any_of(equation.terms, [](const LinearTerm& term) {
              const auto unsupported = term.derivative == Derivative::Z ||
                                       (Dimension == 1 && term.derivative == Derivative::Y);
              return unsupported && term.coefficient != 0.0;
            });
          })) {
        throw std::invalid_argument("Equation derivative exceeds the LinearDG mesh dimension");
      }
    }
  }

  void initialize_flux_matrices() {
    validate_flux_policy();
    const auto fields = equations.fields();
    LOG_DEBUG("Assembling {}D linear DG flux matrices for {} fields", Dimension, fields.size());
    for (auto& matrix : flux_matrices) {
      matrix.assign(fields.size(), std::vector<double>(fields.size(), 0.0));
    }
    constexpr std::array directions = {Derivative::X, Derivative::Y, Derivative::Z};
    for (std::size_t direction = 0; direction < Dimension; ++direction) {
      for (std::size_t equation = 0; equation < fields.size(); ++equation) {
        for (std::size_t field = 0; field < fields.size(); ++field) {
          flux_matrices[direction][equation][field] =
              equations.flux_coefficient(fields[equation], directions[direction], fields[field]);
        }
      }
    }
  }

  void initialize_normal_flux_matrix(FaceFluxOperators& face, std::size_t fields) const {
    face.normal_flux_matrix.resize(fields * fields);
    for (std::size_t equation = 0; equation < fields; ++equation) {
      for (std::size_t field = 0; field < fields; ++field) {
        for (std::size_t direction = 0; direction < Dimension; ++direction) {
          face.normal_flux_matrix[equation * fields + field] +=
              face.normal[direction] * flux_matrices[direction][equation][field];
        }
      }
    }
  }

  void initialize_normal_flux_operators() {
    const auto fields = equations.fields().size();
    LOG_DEBUG("Precomputing DG normal flux operators");
    normal_flux_operators.reserve(integration.cells().size());
    for (const auto& cell_rule : integration.cells()) {
      auto& cell = normal_flux_operators.emplace_back();
      if (!cell_rule.has_group()) continue;

      cell.faces.reserve(cell_rule.faces.size());
      for (const auto& face_rule : cell_rule.faces) {
        auto& face = cell.faces.emplace_back();
        std::copy_n(face_rule.outward_normal.begin(), Dimension, face.normal.begin());
        initialize_normal_flux_matrix(face, fields);
      }
    }
  }

  static void validate_state_strides(const xt::xarray<double>& values) {
    const auto& shape = values.shape();
    const auto& strides = values.strides();
    const auto wrong_field_stride = shape[2] > 1 && strides[2] != 1;
    const auto wrong_node_stride = shape[1] > 1 && std::cmp_not_equal(strides[1], shape[2]);
    const auto wrong_cell_stride =
        shape[0] > 1 && std::cmp_not_equal(strides[0], shape[1] * shape[2]);
    if (wrong_field_stride || wrong_node_stride || wrong_cell_stride) {
      throw std::invalid_argument("DGState arrays must use contiguous row-major storage");
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
      validate_state_strides(values);
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

  static double sample_initial_condition(const InitialCondition& value,
                                         const xt::xarray<double>& nodes, std::size_t cell,
                                         std::size_t node) {
    if constexpr (Dimension == 1) {
      return value(nodes(cell, node, 0));
    } else if constexpr (Dimension == 2) {
      return value(nodes(cell, node, 0), nodes(cell, node, 1));
    } else {
      return value(nodes(cell, node, 0), nodes(cell, node, 1), nodes(cell, node, 2));
    }
  }

  void set_initial_condition(State& state, const Field& field,
                             const InitialCondition& value) const {
    if (!value) throw std::invalid_argument("Initial condition callback cannot be empty");
    validate_state(state);
    const auto field_index = equations.field_index(field);
    const auto groups = space->element_groups();
    LOG_DEBUG("Sampling initial condition for DG field {} across {} groups", field_index,
              groups.size());
    for (std::size_t group_index = 0; group_index < groups.size(); ++group_index) {
      const auto& nodes = groups[group_index].nodes;
      for (std::size_t cell = 0; cell < nodes.shape()[0]; ++cell) {
        for (std::size_t node = 0; node < nodes.shape()[1]; ++node) {
          state.groups[group_index](cell, node, field_index) =
              sample_initial_condition(value, nodes, cell, node);
        }
      }
    }
  }

  State initial_state(const Field& field, const InitialCondition& value) const {
    auto state = zero_state();
    set_initial_condition(state, field, value);
    return state;
  }

  static NodalValues nodal_values(const State& state, const DGIntegration::CellRule& cell,
                                  std::size_t fields) {
    const auto& group = state.groups[cell.group_index];
    const auto node_count = group.shape()[1];
    return {group.data() + cell.group_cell_index * node_count * fields, node_count};
  }

  std::array<double, Dimension> reference_flux_at_node(
      const double* node_values, std::size_t fields, std::size_t equation,
      const std::array<double, 9>& inverse_jacobian) const {
    std::array<double, Dimension> physical_flux{};
    for (std::size_t direction = 0; direction < Dimension; ++direction) {
      for (std::size_t field = 0; field < fields; ++field) {
        physical_flux[direction] += flux_matrices[direction][equation][field] * node_values[field];
      }
    }
    std::array<double, Dimension> transformed_flux{};
    for (std::size_t axis = 0; axis < Dimension; ++axis) {
      for (std::size_t direction = 0; direction < Dimension; ++direction) {
        transformed_flux[axis] += inverse_jacobian[axis * 3 + direction] * physical_flux[direction];
      }
    }
    return transformed_flux;
  }

  static double flux_divergence(const ReferenceFlux& reference_flux,
                                const std::array<const double*, 3>& derivatives,
                                std::size_t node_count, std::size_t fields, std::size_t test,
                                std::size_t equation) {
    double value = 0.0;
    for (std::size_t trial = 0; trial < node_count; ++trial) {
      const auto index = test * node_count + trial;
      for (std::size_t axis = 0; axis < Dimension; ++axis) {
        value += derivatives[axis][index] * reference_flux[axis][trial * fields + equation];
      }
    }
    return value;
  }

  void set_volume_derivative(const DGIntegration::CellRule& cell, NodalValues values,
                             std::size_t fields, double* derivative,
                             ReferenceFlux& reference_flux) const {
    const auto node_count = values.node_count;
    const std::array derivatives = {cell.reference->derivative_r.data(),
                                    cell.reference->derivative_s.data(),
                                    cell.reference->derivative_t.data()};
    for (auto& flux : reference_flux) flux.resize(node_count * fields);
    // The affine inverse Jacobian is constant, so transform the flux before differentiating.
    for (std::size_t trial = 0; trial < node_count; ++trial) {
      for (std::size_t equation = 0; equation < fields; ++equation) {
        const auto flux = reference_flux_at_node(values.data + trial * fields, fields, equation,
                                                 cell.inverse_jacobian);
        for (std::size_t axis = 0; axis < Dimension; ++axis) {
          reference_flux[axis][trial * fields + equation] = flux[axis];
        }
      }
    }
    for (std::size_t test = 0; test < node_count; ++test) {
      for (std::size_t equation = 0; equation < fields; ++equation) {
        derivative[test * fields + equation] =
            -flux_divergence(reference_flux, derivatives, node_count, fields, test, equation);
      }
    }
  }

  static void interpolate_trace(std::vector<double>& trace, const double* basis,
                                NodalValues values) {
    const auto fields = trace.size();
    for (std::size_t node = 0; node < values.node_count; ++node) {
      for (std::size_t field = 0; field < fields; ++field) {
        trace[field] += basis[node] * values.data[node * fields + field];
      }
    }
  }

  std::vector<double> face_flux_correction(const FaceFluxOperators& operators,
                                           FaceScratch& scratch) const {
    const auto fields = scratch.left.size();
    for (std::size_t equation = 0; equation < fields; ++equation) {
      const auto* normal_matrix = operators.normal_flux_matrix.data() + equation * fields;
      double interior_flux = 0.0;
      double exterior_flux = 0.0;
      for (std::size_t field = 0; field < fields; ++field) {
        interior_flux += normal_matrix[field] * scratch.left[field];
        exterior_flux += normal_matrix[field] * scratch.right[field];
      }
      scratch.left_flux[equation] = interior_flux;
      scratch.right_flux[equation] = exterior_flux;
    }
    auto correction = numerical_flux(scratch.left, scratch.right, operators.normal,
                                     scratch.left_flux, scratch.right_flux);
    if (correction.size() != fields) {
      throw std::invalid_argument("Numerical flux must return one value per equation field");
    }
    for (std::size_t equation = 0; equation < fields; ++equation) {
      correction[equation] -= scratch.left_flux[equation];
    }
    return correction;
  }

  void add_face_derivative(const DGIntegration::FaceRule& face, const FaceFluxOperators& operators,
                           NodalValues cell, NodalValues neighbor, double* derivative,
                           FaceScratch& scratch) const {
    const auto fields = scratch.left.size();
    const auto node_count = cell.node_count;
    const auto face_points = face.point_count();
    for (std::size_t q = 0; q < face_points; ++q) {
      std::ranges::fill(scratch.left, 0.0);
      interpolate_trace(scratch.left, face.minus->trace.data() + q * node_count, cell);
      if (face.is_boundary()) {
        std::ranges::copy(scratch.left, scratch.right.begin());
        if (boundary_state) boundary_state(scratch.left, operators.normal, scratch.right);
      } else {
        std::ranges::fill(scratch.right, 0.0);
        interpolate_trace(scratch.right, face.plus->trace.data() + q * neighbor.node_count,
                          neighbor);
      }

      const auto correction = face_flux_correction(operators, scratch);
      const auto* lift = face.minus->lift.data() + q * node_count;
      for (std::size_t test = 0; test < node_count; ++test) {
        for (std::size_t equation = 0; equation < fields; ++equation) {
          derivative[test * fields + equation] -=
              face.lift_scale * lift[test] * correction[equation];
        }
      }
    }
  }

  State rhs(const State& state) const {
    validate_state(state);
    auto output = zero_state();
    const auto fields = equations.fields().size();
    FaceScratch scratch(fields);
    ReferenceFlux reference_flux;

    const auto cell_rules = integration.cells();
    for (std::size_t cell_index = 0; cell_index < cell_rules.size(); ++cell_index) {
      const auto& cell = cell_rules[cell_index];
      if (!cell.has_group()) continue;
      const auto values = nodal_values(state, cell, fields);
      auto* derivative = output.groups[cell.group_index].data() +
                         cell.group_cell_index * values.node_count * fields;
      set_volume_derivative(cell, values, fields, derivative, reference_flux);

      const auto& operators = normal_flux_operators[cell_index];
      for (std::size_t face_index = 0; face_index < cell.faces.size(); ++face_index) {
        const auto& face = cell.faces[face_index];
        const auto& neighbor = face.is_boundary() ? cell : integration.cell(face.neighbor_cell);
        add_face_derivative(face, operators.faces[face_index], values,
                            nodal_values(state, neighbor, fields), derivative, scratch);
      }
    }
    return output;
  }
};

template <std::size_t Dimension>
LinearDG<Dimension>::LinearDG(mesh::Mesh& mesh, const std::vector<unsigned>& orders,
                              EquationSystem equations, WaveSpeed wave_speed,
                              BoundaryState boundary_state)
    : LinearDG(mesh, orders, std::move(equations), rusanov(std::move(wave_speed)),
               std::move(boundary_state)) {}

template <std::size_t Dimension>
LinearDG<Dimension>::LinearDG(mesh::Mesh& mesh, const std::vector<unsigned>& orders,
                              EquationSystem equations, NumericalFlux numerical_flux,
                              BoundaryState boundary_state)
    : m_impl(std::make_unique<Impl>(mesh, orders, std::move(equations), std::move(numerical_flux),
                                    std::move(boundary_state))) {}

template <std::size_t Dimension>
LinearDG<Dimension>::LinearDG(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
                              WaveSpeed wave_speed, BoundaryState boundary_state)
    : LinearDG(mesh, space, std::move(equations), rusanov(std::move(wave_speed)),
               std::move(boundary_state)) {}

template <std::size_t Dimension>
LinearDG<Dimension>::LinearDG(mesh::Mesh& mesh, const DGSpace& space, EquationSystem equations,
                              NumericalFlux numerical_flux, BoundaryState boundary_state)
    : m_impl(std::make_unique<Impl>(mesh, space, std::move(equations), std::move(numerical_flux),
                                    std::move(boundary_state))) {}

template <std::size_t Dimension>
LinearDG<Dimension>::~LinearDG() = default;

template <std::size_t Dimension>
LinearDG<Dimension>::LinearDG(LinearDG&&) noexcept = default;

template <std::size_t Dimension>
LinearDG<Dimension>& LinearDG<Dimension>::operator=(LinearDG&&) noexcept = default;

template <std::size_t Dimension>
const DGSpace& LinearDG<Dimension>::space() const {
  return *m_impl->space;
}

template <std::size_t Dimension>
const EquationSystem& LinearDG<Dimension>::equations() const {
  return m_impl->equations;
}

template <std::size_t Dimension>
auto LinearDG<Dimension>::rusanov(WaveSpeed wave_speed) -> NumericalFlux {
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

template <std::size_t Dimension>
auto LinearDG<Dimension>::zero_state() const -> State {
  return m_impl->zero_state();
}

template <std::size_t Dimension>
auto LinearDG<Dimension>::initial_state(const Field& field, const InitialCondition& value) const
    -> State {
  return m_impl->initial_state(field, value);
}

template <std::size_t Dimension>
void LinearDG<Dimension>::set_initial_condition(State& state, const Field& field,
                                                const InitialCondition& value) const {
  m_impl->set_initial_condition(state, field, value);
}

template <std::size_t Dimension>
auto LinearDG<Dimension>::rhs(const State& state) const -> State {
  return m_impl->rhs(state);
}

template class LinearDG<1>;
template class LinearDG<2>;
template class LinearDG<3>;

}  // namespace oiseau::dg::solver
