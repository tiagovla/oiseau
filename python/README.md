# Oiseau Python bindings

From `oiseau/python/`, sync the environment with CPython 3.13 or later and build the bindings.
The build requires Git, CMake 3.26 or later, a C++23 compiler, and BLAS/LAPACK development
libraries. CMake fetches the C++ dependencies on the first build.

```sh
uv sync --locked
uv run --locked python -c 'from oiseau.dg import LinearDG2D; print(LinearDG2D.__name__)'
```

Use `uv sync --locked --extra gmsh` to install the optional Gmsh bindings.

Run the 2D perfectly conducting cavity example from the same directory to plot
the electric field after half a period.

```sh
uv run python -m demo.dg_pec_cavity_2d
```

To save the plot without a display, pass `--output cavity_2d.png`.

After editing the C++ bindings, rebuild with `uv sync --reinstall-package oiseau`.
