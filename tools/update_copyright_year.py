#!/usr/bin/env python3
# Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Update copyright years in files passed by pre-commit."""

import datetime
import re
import sys
from pathlib import Path


CURRENT_YEAR = datetime.date.today().year
COPYRIGHT = re.compile(
    r"(Copyright \(C\) )(\d{4})(?:-(\d{4}))?( Tiago V\. L\. Amorim \(@tiagovla\))"
)


def update_copyright_year(content: str) -> str:
    def replace(match: re.Match[str]) -> str:
        return f"{match.group(1)}{CURRENT_YEAR}{match.group(4)}"

    return COPYRIGHT.sub(replace, content)


def main() -> None:
    for filename in sys.argv[1:]:
        path = Path(filename)
        try:
            content = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue

        updated_content = update_copyright_year(content)
        if updated_content != content:
            path.write_text(updated_content, encoding="utf-8")


if __name__ == "__main__":
    main()
