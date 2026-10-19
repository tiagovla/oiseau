# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys
import subprocess
from pathlib import Path

import ninja
from setuptools import Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = ""):
        super().__init__(name, sources=[])
        self.sourcedir = str(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        build_temp = Path(self.build_temp) / ext.name
        build_temp.mkdir(parents=True, exist_ok=True)

        cmake_command = [
            "cmake",
            ext.sourcedir,
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DCMAKE_BUILD_TYPE=Release",
            "-GNinja",
            f"-DCMAKE_MAKE_PROGRAM={Path(ninja.BIN_DIR) / 'ninja'}",
        ]

        if "CMAKE_ARGS" in os.environ:
            cmake_command += os.environ["CMAKE_ARGS"].split()

        cmake_build_command = ["cmake", "--build", ".", "-j8"]

        try:
            subprocess.check_call(cmake_command, cwd=build_temp)
            subprocess.check_call(cmake_build_command, cwd=build_temp)
        except subprocess.CalledProcessError as e:
            print(f"\n❌ Error while building {ext.name}: {e}\n", file=sys.stderr)
            raise


def build(setup_kwargs):
    os.environ["CFLAGS"] = "-O3"

    setup_kwargs.update(
        {
            "version": "0.1.0",
            "author": "tiagovla",
            "author_email": "tiagovla@gmail.com",
            "description": "DG Framework",
            "packages": ["oiseau", "oiseau.wrappers"],
            "package_data": {"oiseau.cpp": ["*.h"]},
            "ext_modules": [CMakeExtension("oiseau.cpp")],
            "cmdclass": {"build_ext": CMakeBuild},
            "zip_safe": False,
            "setup_requires": ["pybind11"],
        }
    )
