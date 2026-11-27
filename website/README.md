<!--
 Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)

 This file is part of oiseau (https://github.com/tiagovla/oiseau)

 SPDX-License-Identifier: GPL-3.0-or-later
-->

# Oiseau website

Install Node.js 22.23.4 or later, then run these commands from the repository root.

```sh
npm ci --prefix website
npm run dev --prefix website
```

Astro serves the landing page at `http://localhost:4321/`. Run `make serve` to
build and serve the full site with the C++ API documentation and any enabled
coverage report. The CMake `website` target writes the combined site to
`build/web/` with the default build directory.

The GitHub Pages workflow sets `OISEAU_SITE_BASE=/oiseau` when building the site.

## Theory topics

Add an Astro file with the topic text in `src/content/theory/`, then register it
in `src/content/theory/topics.ts`. The page uses the list for its section order,
headings, and table of contents. Use a unique `id` for each section link.
