// Copyright (C) 2026 Tiago V. L. Amorim (@tiagovla)
//
// This file is part of oiseau (https://github.com/tiagovla/oiseau)
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <benchmark/benchmark.h>

#include <cstddef>
#include <random>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xmath.hpp>
#include <xtensor/core/xtensor_forward.hpp>
#include <xtensor/generators/xbuilder.hpp>
#include <xtensor/generators/xrandom.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/views/xview.hpp>

using std::size_t;

constexpr std::size_t SIZE = 1000;

xt::xarray<double> generate_xarray(std::size_t seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  xt::xarray<double> a = xt::zeros<double>({SIZE, SIZE});
  for (std::size_t i = 0; i < SIZE; ++i) {
    for (std::size_t j = 0; j < SIZE; ++j) {
      a(i, j) = dist(gen);
    }
  }
  return a;
}

xt::xtensor<double, 2> generate_xtensor(std::size_t seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  xt::xtensor<double, 2> a = xt::zeros<double>({SIZE, SIZE});
  for (std::size_t i = 0; i < SIZE; ++i) {
    for (size_t j = 0; j < SIZE; ++j) {
      a(i, j) = dist(gen);
    }
  }
  return a;
}

// --------------------- Allocation ---------------------
static void BM_Xarray_Allocation(benchmark::State& state) {
  for (auto _ : state) {
    xt::xarray<double> a = xt::zeros<double>({SIZE, SIZE});
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_Xarray_Allocation)->Iterations(1000);

static void BM_Xtensor_Allocation(benchmark::State& state) {
  for (auto _ : state) {
    xt::xtensor<double, 2> a = xt::zeros<double>({SIZE, SIZE});
    benchmark::DoNotOptimize(a);
  }
}
BENCHMARK(BM_Xtensor_Allocation)->Iterations(1000);

// --------------------- Addition ---------------------
static void BM_Xarray_Addition(benchmark::State& state) {
  auto a = generate_xarray(42);
  auto b = generate_xarray(43);
  for (auto _ : state) {
    xt::xarray<double> c = a + b;
    benchmark::DoNotOptimize(c);
  }
}
BENCHMARK(BM_Xarray_Addition)->Iterations(100);

static void BM_Xtensor_Addition(benchmark::State& state) {
  auto a = generate_xtensor(42);
  auto b = generate_xtensor(43);
  for (auto _ : state) {
    xt::xtensor<double, 2> c = a + b;
    benchmark::DoNotOptimize(c);
  }
}
BENCHMARK(BM_Xtensor_Addition)->Iterations(100);

// --------------------- Sum ---------------------
static void BM_Xarray_Sum(benchmark::State& state) {
  auto a = generate_xarray(44);
  for (auto _ : state) {
    double s = xt::sum(a)();
    benchmark::DoNotOptimize(s);
  }
}
BENCHMARK(BM_Xarray_Sum)->Iterations(1000);

static void BM_Xtensor_Sum(benchmark::State& state) {
  auto a = generate_xtensor(44);
  for (auto _ : state) {
    double s = xt::sum(a)();
    benchmark::DoNotOptimize(s);
  }
}
BENCHMARK(BM_Xtensor_Sum)->Iterations(1000);

// --------------------- Slicing ---------------------
static void BM_Xarray_Slicing(benchmark::State& state) {
  auto a = generate_xarray(45);
  for (auto _ : state) {
    auto view = xt::view(a, xt::range(100, 200), xt::range(100, 200));
    benchmark::DoNotOptimize(view);
  }
}
BENCHMARK(BM_Xarray_Slicing)->Iterations(1000000);

static void BM_Xtensor_Slicing(benchmark::State& state) {
  auto a = generate_xtensor(45);
  for (auto _ : state) {
    auto view = xt::view(a, xt::range(100, 200), xt::range(100, 200));
    benchmark::DoNotOptimize(view);
  }
}
BENCHMARK(BM_Xtensor_Slicing)->Iterations(1000000);

BENCHMARK_MAIN();
