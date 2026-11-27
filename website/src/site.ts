// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

export const ownerUrl = 'https://github.com/tiagovla';
export const repositoryUrl = `${ownerUrl}/oiseau`;

const base = import.meta.env.BASE_URL.replace(/\/$/, '');

export function sitePath(path: string): string {
  return `${base}/${path}`;
}
