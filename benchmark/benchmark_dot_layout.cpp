// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <benchmark/benchmark.h>
#include <omp.h>

// #undef WITH_OPENBLAS to disable openblas

#include <cstddef>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xlayout.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xrandom.hpp>
#include <xtensor/io/xio.hpp>

#define RANGE 2 << 8, 2 << 20
#define MULTIPLIER 16
constexpr std::size_t N = 32;

// Benchmark: A row-major, B row-major
void Dot_Arow_Brow(benchmark::State& state) {
  using namespace xt;
  xtensor<double, 2> a = random::rand<double>({N, N});
  std::vector<std::size_t> shape = {N, static_cast<std::size_t>(state.range(0))};
  xtensor<double, 2> b = random::rand<double>(shape);

  for (auto _ : state) {
    xtensor<double, 2> res = xt::linalg::dot(a, b);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(Dot_Arow_Brow)->RangeMultiplier(MULTIPLIER)->Range(RANGE)->Unit(benchmark::kMillisecond);

// Benchmark: A row-major, B col-major
void Dot_Arow_Bcol(benchmark::State& state) {
  using namespace xt;

  xtensor<double, 2, layout_type::row_major> a = random::rand<double>({N, N});
  std::vector<std::size_t> shape = {N, static_cast<std::size_t>(state.range(0))};
  xtensor<double, 2, layout_type::column_major> b = random::rand<double>(shape);

  for (auto _ : state) {
    xtensor<double, 2> res = xt::linalg::dot(a, b);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(Dot_Arow_Bcol)->RangeMultiplier(MULTIPLIER)->Range(RANGE)->Unit(benchmark::kMillisecond);

// Benchmark: A col-major, B row-major
void Dot_Acol_Brow(benchmark::State& state) {
  using namespace xt;

  xtensor<double, 2, layout_type::column_major> a = random::rand<double>({N, N});
  std::vector<std::size_t> shape = {N, static_cast<std::size_t>(state.range(0))};
  xtensor<double, 2, layout_type::row_major> b = random::rand<double>(shape);

  for (auto _ : state) {
    xtensor<double, 2> res = xt::linalg::dot(a, b);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(Dot_Acol_Brow)->RangeMultiplier(MULTIPLIER)->Range(RANGE)->Unit(benchmark::kMillisecond);

// Benchmark: A col-major, B col-major
void Dot_Acol_Bcol(benchmark::State& state) {
  using namespace xt;

  xtensor<double, 2, layout_type::column_major> a = random::rand<double>({N, N});
  std::vector<std::size_t> shape = {N, static_cast<std::size_t>(state.range(0))};
  xtensor<double, 2, layout_type::column_major> b = random::rand<double>(shape);

  for (auto _ : state) {
    xtensor<double, 2> res = xt::linalg::dot(a, b);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(Dot_Acol_Bcol)->RangeMultiplier(MULTIPLIER)->Range(RANGE)->Unit(benchmark::kMillisecond);

int main(int argc, char** argv) {
  omp_set_num_threads(8);  // Set the number of threads to 8
  benchmark::MaybeReenterWithoutASLR(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
