#!/usr/bin/env python3
# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Extend copyright year ranges in staged files passed by pre-commit."""

import datetime
import os
import re
import subprocess
import sys
from pathlib import Path

CURRENT_YEAR = datetime.date.today().year
COPYRIGHT = re.compile(
    r"(Copyright \(C\) )(\d{4})(?:-(\d{4}))?( Tiago V\. L\. Amorim \(@tiagovla\))"
)


def update_copyright_year(content: str) -> str:
    def replace(match: re.Match[str]) -> str:
        first_year = match.group(2)
        last_year = int(match.group(3) or first_year)
        if CURRENT_YEAR <= last_year:
            return match.group(0)
        return f"{match.group(1)}{first_year}-{CURRENT_YEAR}{match.group(4)}"

    return COPYRIGHT.sub(replace, content)


def main() -> None:
    if len(sys.argv) == 1:
        return

    output = subprocess.check_output(
        [
            "git",
            "diff",
            "--cached",
            "--relative",
            "--name-only",
            "--diff-filter=ACMRT",
            "-z",
        ]
    )
    staged_paths = {Path(os.fsdecode(name)).absolute() for name in output.split(b"\0") if name}
    for filename in sys.argv[1:]:
        path = Path(filename).absolute()
        if path not in staged_paths or path.is_symlink():
            continue
        try:
            content = path.read_text(encoding="utf-8")
        except (FileNotFoundError, UnicodeDecodeError):
            continue

        updated_content = update_copyright_year(content)
        if updated_content != content:
            path.write_text(updated_content, encoding="utf-8")


if __name__ == "__main__":
    main()
