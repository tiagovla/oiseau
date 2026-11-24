# Copyright (C) 2025 Tiago V. L. Amorim (@tiagovla)
#
# This file is part of oiseau (https://github.com/tiagovla/oiseau)
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Use Go image to build the binary
FROM golang:1.27 AS build

# Build the same addlicense release used by pre-commit.
RUN CGO_ENABLED=0 GOOS=linux GOBIN=/app go install github.com/google/addlicense@v1.2.0

# Minimal final image
FROM scratch

# Source to be scanned should be mounted to /src
WORKDIR /src

# Copy only the built binary from builder
COPY --from=build /app/addlicense /app/addlicense

ENTRYPOINT ["/app/addlicense"]
