// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

import ConservationLaws from './ConservationLaws.astro';
import ReferenceElements from './ReferenceElements.astro';

export const topics = [
  {
    id: 'conservation-laws',
    title: 'Conservation laws',
    Content: ConservationLaws,
  },
  {
    id: 'reference-elements',
    title: 'Reference elements',
    Content: ReferenceElements,
  },
];
