# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Sao_Paulo

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang-tidy \
    cmake \
    cppcheck \
    curl \
    doxygen \
    gcovr \
    gfortran \
    git \
    gnupg \
    liblapack-dev \
    libomp-dev \
    libopenblas-dev \
    ninja-build \
    pkg-config \
    python3 \
    software-properties-common \
    unzip \
    zip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
