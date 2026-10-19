# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import gmsh

gmsh.initialize()
gmsh.model.add("split_triangles_quads")

lc = 0.5
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

s_left = gmsh.model.occ.addRectangle(0, 0, 0, 0.5, 1)
s_right = gmsh.model.occ.addRectangle(0.5, 0, 0, 0.5, 1)

gmsh.model.occ.synchronize()

pg_triangles = gmsh.model.addPhysicalGroup(2, [s_left], name="Region_Triangles")
pg_quads = gmsh.model.addPhysicalGroup(2, [s_right], name="Region_Quads")

gmsh.model.mesh.setRecombine(2, s_right, True)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
gmsh.option.setNumber("Mesh.Algorithm", 6)  # 2D Frontal-Delaunay

gmsh.model.occ.removeAllDuplicates()

gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
gmsh.write("mesh2d_split_triangles_quads.msh")

# gmsh.fltk.run()

gmsh.finalize()
