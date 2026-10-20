// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_geometry.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/plotting/triplot.hpp"

namespace {

void print_matrix(const std::array<double, 9>& matrix, std::size_t dimension) {
  const std::array<std::size_t, 2> shape{dimension, dimension};
  const std::array<std::size_t, 2> strides{3, 1};
  std::cout << xt::adapt(matrix, shape, strides);
}

void plot_normals(oiseau::mesh::Mesh& mesh, const oiseau::dg::DGGeometry& geometry) {
  const auto cell_types = mesh.topology().cell_types();
  const auto connectivity = mesh.topology().conn();
  const auto coordinates = mesh.geometry().x();
  const auto coordinate_dimension = mesh.geometry().dim();
  const auto dimension = static_cast<std::size_t>(mesh.topology().dimension());
  const std::vector<std::size_t> coordinate_shape = {coordinates.size() / coordinate_dimension,
                                                     coordinate_dimension};
  const auto coordinate_matrix =
      xt::adapt(coordinates.data(), coordinates.size(), xt::no_ownership(), coordinate_shape);

  std::size_t face_count = 0;
  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    if (geometry.has_cell(cell_index)) face_count += geometry.cell(cell_index).faces.size();
  }
  xt::xarray<double> face_centers = xt::zeros<double>({face_count, dimension});
  xt::xarray<double> face_normals = xt::zeros<double>({face_count, dimension});
  std::size_t face_row = 0;
  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    if (!geometry.has_cell(cell_index)) continue;
    const auto& cell_geometry = geometry.cell(cell_index);
    const auto local_faces =
        cell_types[cell_index]->get_entity_vertices(cell_geometry.dimension - 1);
    for (std::size_t face_index = 0; face_index < local_faces.size(); ++face_index) {
      std::vector<std::size_t> face_vertices;
      face_vertices.reserve(local_faces[face_index].size());
      for (const auto local_vertex : local_faces[face_index]) {
        face_vertices.push_back(connectivity[cell_index][static_cast<std::size_t>(local_vertex)]);
      }
      const auto physical_face = xt::view(coordinate_matrix, xt::keep(face_vertices), xt::all());
      const auto center = xt::mean(physical_face, {0});
      for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
        face_centers(face_row, coordinate) = center(coordinate);
        face_normals(face_row, coordinate) =
            cell_geometry.faces[face_index].outward_normal[coordinate];
      }
      ++face_row;
    }
  }

  plt::scoped_interpreter interpreter{};
  if (mesh.topology().dimension() == 2) {
    auto [figure, axes] =
        plt::subplots(1, 1, "figsize"_a = py::make_tuple(7.0, 6.0), "constrained_layout"_a = true);
    oiseau::plotting::triplot(axes, mesh);
    axes.quiver(xt::col(face_centers, 0), xt::col(face_centers, 1), xt::col(face_normals, 0),
                xt::col(face_normals, 1), "angles"_a = "xy", "scale_units"_a = "xy",
                "scale"_a = 8.0, "color"_a = "tab:red", "width"_a = 0.004);
    axes.set_aspect("equal");
    axes.set_xlabel("x");
    axes.set_ylabel("y");
    axes.set_title("DGGeometry outward face normals");
  } else {
    auto figure = plt::figure("figsize"_a = py::make_tuple(7.0, 6.0));
    auto axes = figure.add_subplot(111, "projection"_a = "3d");
    for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
      if (!geometry.has_cell(cell_index)) continue;
      const auto local_faces = cell_types[cell_index]->get_entity_vertices(2);
      for (const auto& local_face : local_faces) {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        for (const auto local_vertex : local_face) {
          const auto vertex = connectivity[cell_index][static_cast<std::size_t>(local_vertex)];
          x.push_back(coordinates[vertex * coordinate_dimension]);
          y.push_back(coordinates[vertex * coordinate_dimension + 1]);
          z.push_back(coordinates[vertex * coordinate_dimension + 2]);
        }
        x.push_back(x.front());
        y.push_back(y.front());
        z.push_back(z.front());
        axes.plot(x, y, z, "color"_a = "black");
      }
    }
    axes.quiver(xt::col(face_centers, 0), xt::col(face_centers, 1), xt::col(face_centers, 2),
                xt::col(face_normals, 0), xt::col(face_normals, 1), xt::col(face_normals, 2),
                "length"_a = 0.15, "normalize"_a = true, "color"_a = "tab:red");
    axes.set_xlabel("x");
    axes.set_ylabel("y");
    axes.set_zlabel("z");
    axes.set_title("DGGeometry outward face normals");
  }
  plt::show();
}

}  // namespace

int main(int argc, char** argv) {
  using namespace oiseau;

  const std::filesystem::path mesh_path =
      argc > 1 ? argv[1] : "demo/meshes/mesh2d_two_triangles.msh";
  auto mesh = io::gmsh_read_from_path(mesh_path);
  const dg::DGGeometry geometry(mesh);

  std::cout << "DGGeometry demo\n"
            << "Mesh: " << mesh_path << "\n"
            << "Each affine cell stores J = d(x)/d(r), J^-1, and outward face data.\n\n";

  for (std::size_t cell_index = 0; cell_index < mesh.topology().n_cells(); ++cell_index) {
    if (!geometry.has_cell(cell_index)) continue;
    const auto& cell = geometry.cell(cell_index);
    std::cout << "Mesh cell " << cell.mesh_cell_index << "\n  J     = \n";
    print_matrix(cell.jacobian, cell.dimension);
    std::cout << "\n  J^-1  = ";
    print_matrix(cell.inverse_jacobian, cell.dimension);
    std::cout << "\n  det J = " << cell.determinant << '\n';
    for (std::size_t face_index = 0; face_index < cell.faces.size(); ++face_index) {
      const auto& face = cell.faces[face_index];
      std::cout << "  face " << face_index << ": normal = (";
      for (std::size_t coordinate = 0; coordinate < cell.dimension; ++coordinate) {
        if (coordinate != 0) std::cout << ", ";
        std::cout << face.outward_normal[coordinate];
      }
      std::cout << "), measure = " << face.measure
                << ", surface Jacobian = " << face.surface_jacobian << '\n';
    }
    std::cout << '\n';
  }

  plot_normals(mesh, geometry);

  return 0;
}
