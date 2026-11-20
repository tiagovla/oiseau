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
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = ""):
        super().__init__(name, sources=[], py_limited_api=True)
        self.sourcedir = str(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    def get_ext_filename(self, fullname: str) -> str:
        filename = super().get_ext_filename(fullname)
        stem, marker, platform_suffix = filename.partition(".abi3-")
        if marker:
            # Nanobind uses the unqualified suffix understood by Python 3.13 and 3.14.
            return f"{stem}.abi3{Path(platform_suffix).suffix}"
        return filename

    def build_extension(self, ext):
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        build_temp = Path(self.build_temp) / ext.name
        build_temp.mkdir(parents=True, exist_ok=True)

        cmake_command = [
            "cmake",
            ext.sourcedir,
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPython_EXECUTABLE={sys.executable}",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
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


setup(
    ext_modules=[CMakeExtension("oiseau.cpp")],
    cmdclass={"build_ext": CMakeBuild},
    options={"bdist_wheel": {"py_limited_api": "cp313"}},
)
