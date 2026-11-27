// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

import { defineConfig } from 'astro/config';
import tailwindcss from 'tailwindcss';

export default defineConfig({
  site: 'https://tiagovla.github.io',
  base: process.env.OISEAU_SITE_BASE ?? '/',
  vite: {
    css: {
      postcss: {
        plugins: [
          tailwindcss({
            content: ['./src/**/*.astro'],
            theme: {
              extend: {
                colors: {
                  'oiseau-dark-bg': '#0F2027',
                  'oiseau-deep-blue': '#0A192F',
                  'oiseau-light-text': '#CCD6F6',
                  'oiseau-teal-accent': '#64FFDA',
                  'oiseau-orange-accent': '#FF8C00',
                  'oiseau-nav-link': '#8892b0',
                },
                fontFamily: {
                  sans: [
                    'Inter',
                    'system-ui',
                    '-apple-system',
                    'BlinkMacSystemFont',
                    'Segoe UI',
                    'Roboto',
                    'Helvetica Neue',
                    'Arial',
                    'Noto Sans',
                    'sans-serif',
                    'Apple Color Emoji',
                    'Segoe UI Emoji',
                    'Segoe UI Symbol',
                    'Noto Color Emoji',
                  ],
                },
              },
            },
          }),
        ],
      },
    },
  },
});
