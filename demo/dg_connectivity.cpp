// Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "oiseau/dg/dg_connectivity.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "oiseau/io/gmsh.hpp"
#include "oiseau/mesh/cell.hpp"
#include "oiseau/mesh/mesh.hpp"

namespace {

void print_vertex(std::size_t vertex, const oiseau::mesh::Geometry& geometry) {
  const auto coordinates = geometry.x();
  std::cout << vertex << " (";
  for (std::size_t coordinate = 0; coordinate < geometry.dim(); ++coordinate) {
    if (coordinate != 0) std::cout << ", ";
    std::cout << coordinates[vertex * geometry.dim() + coordinate];
  }
  std::cout << ')';
}

void print_face_vertices(const oiseau::mesh::Cell& cell,
                         const std::vector<std::size_t>& cell_vertices, std::size_t face_index,
                         const oiseau::mesh::Geometry& geometry) {
  const auto local_faces = cell.get_entity_vertices(cell.dimension() - 1);
  std::cout << "    vertices: ";
  for (std::size_t vertex = 0; vertex < local_faces[face_index].size(); ++vertex) {
    if (vertex != 0) std::cout << ", ";
    const auto local_vertex = local_faces[face_index][vertex];
    if (local_vertex < 0 || static_cast<std::size_t>(local_vertex) >= cell_vertices.size()) {
      throw std::invalid_argument("Cell connectivity does not match its cell type");
    }
    print_vertex(cell_vertices[static_cast<std::size_t>(local_vertex)], geometry);
  }
  std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  using namespace oiseau;

  const std::filesystem::path mesh_path =
      argc > 1 ? argv[1] : "demo/meshes/mesh2d_two_triangles_one_quad.msh";
  auto mesh = io::gmsh_read_from_path(mesh_path);
  const auto cell_types = mesh.topology().cell_types();
  const auto mesh_dimension = mesh.topology().dimension();
  if (mesh_dimension != 2) {
    throw std::invalid_argument("This DGConnectivity demo requires a 2D mesh");
  }

  // DGConnectivity uses Topology to find neighboring cells and then adds the
  // vertex permutation needed to align values on an interior face.
  const dg::DGConnectivity connectivity(mesh);
  const auto cell_vertices = mesh.topology().conn();
  const auto& geometry = mesh.geometry();

  std::cout << "DGConnectivity demo\n"
            << "Mesh: " << mesh_path << "\n"
            << "A boundary face has no neighbor. An interior face identifies its\n"
            << "neighbor and the ordering needed to map its vertices to that face.\n\n";

  for (std::size_t cell_index = 0; cell_index < cell_types.size(); ++cell_index) {
    const auto* cell = cell_types[cell_index];
    if (cell->dimension() != mesh_dimension) continue;

    std::cout << "Mesh cell " << cell_index << " (" << cell->name() << ")\n";
    const auto faces = connectivity.faces(cell_index);
    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
      std::cout << "  local face " << face_index << "\n";
      print_face_vertices(*cell, cell_vertices[cell_index], face_index, geometry);

      const auto& connection = faces[face_index];
      if (connection.is_boundary()) {
        std::cout << "    connection: boundary\n";
        continue;
      }

      std::cout << "    connection: mesh cell " << connection.neighbor_cell << ", local face "
                << connection.neighbor_face << "\n"
                << "    vertex permutation: [";
      for (std::size_t vertex = 0; vertex < connection.vertex_permutation.size(); ++vertex) {
        if (vertex != 0) std::cout << ", ";
        std::cout << connection.vertex_permutation[vertex];
      }
      std::cout << "]\n";
    }
    std::cout << '\n';
  }

  return 0;
}
