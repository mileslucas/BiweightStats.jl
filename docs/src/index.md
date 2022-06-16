```@meta
CurrentModule = BiweightStats
```

# BiweightStats.jl

[![Code](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/mileslucas/BiweightStats.jl)
[![Build Status](https://github.com/mileslucas/BiweightStats.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mileslucas/BiweightStats.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BiweightStats.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/mileslucas/BiweightStats.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mileslucas/BiweightStats.jl)
[![License](https://img.shields.io/github/license/mileslucas/BiweightStats.jl?color=yellow)](https://github.com/mileslucas/BiweightStats.jl/blob/main/LICENSE)

Robust statistics based on the biweight transform.

## Installation

BiweightStats.jl is a registered package and can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

pkg> add BiweightStats
```


To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using BiweightStats
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with
For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Usage

The following examples show the biweight [`location`](@ref) and [`scale`](@ref) for a few distributions, and comparing them to the mean and standard deviation.

To run these examples, make sure the following packages need to be installed

```julia
using Pkg; Pkg.add(["Distributions", "StatsPlots"])
```

```@example dists
using BiweightStats
using Distributions: Normal, Logistic, Cauchy, Laplace
using Statistics
using StatsPlots
using Random
rng = Random.seed!(1001)
```

for all our examples, we'll sample 10,000 points

```@example dists
N = 10000
```

```@example dists
function sample_measure_and_plot_dists(rng, dist, N; kwargs...)
    samples = rand(rng, dist, N)

    mu = mean(samples)
    sig = std(samples; mean=mu)

    loc = location(samples; c=9)
    sca = scale(samples; c=9, M=loc)

    p = histogram(samples; normalize=true, fill=0.3, c=1, lab="", leg=:topleft, kwargs...)
    plot!(dist; c=1, lab="", lw=3)
    vline!([mu loc]; c=[2 3], lab=["mean ± std" "biweight location ± scale"], lw=2)
    vline!([mu loc] .+ [-sig -sca; sig sca]; c=[2 3], ls=:dash, lab="", lw=2)
    return p
end
```

### Gaussian

```@example dists
dist = Normal()
sample_measure_and_plot_dists(rng, dist, N; title="Gaussian")
```

### Logistic

```@example dists
dist = Logistic()
sample_measure_and_plot_dists(rng, dist, N; title="Logistic")
```

### Cauchy

```@example dists
dist = Cauchy()
sample_measure_and_plot_dists(rng, dist, N; title="Cauchy", xlim=(-40, 40))
```

### Laplace

```@example dists
dist = Laplace()
sample_measure_and_plot_dists(rng, dist, N; title="Laplace")
```

## Benchmarks

This package has been benchmarked against [astropy.stats](https://github.com/astropy/astropy). The benchmarking code can be found in [`bench/`](https://github.com/mileslucas/BiweightStats.jl/tree/main/bench).

**System Information**

Benchmarks were ran on a 2021 M1 Macbook Pro

```
Julia Version 1.8.0-beta2
Commit b655b2c008 (2022-03-21 12:50 UTC)
Platform Info:
  OS: macOS (arm64-apple-darwin21.3.0)
  CPU: 10 × Apple M1 Max
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, apple-m1)
  Threads: 8 on 8 virtual cores
```

The benchmarks that were ran generated 2 sets of ``n`` normally-distributed samples without any outliers or non-finite numbers added. One of these sets was used for testing both this package's and astropy's implementations of the univariate statistics. For the midcovariance statistic, the data was combined into an ``(N, 2)`` matrix to match the astropy function signature.

```@example bench
using BiweightStats, CSV, DataFrames, StatsPlots # hide
benchdir = joinpath(dirname(pathof(BiweightStats)), "..", "bench") # hide
results = CSV.read(joinpath(benchdir, "benchmark_results.csv"), DataFrame) # hide
shapes = [:o, :square, :diamond, :utriangle, :dtriangle] # hide
groups = groupby(results, :name) # hide
p = plot(xscale=:log10, yscale=:log10, xlabel="n", ylabel="t [s]", # hide
    leg=:topleft, legendfontsize=7) # hide
for (s, (k, g)) in zip(shapes, pairs(groups)) # hide
    @df g scatter!(:n, :t_jl, lab=k.name * " (BiweightStats.jl)", marker=s, c=1) # hide
    @df g scatter!(:n, :t_py, lab=k.name * " (astropy)", marker=s, c=2) # hide
end # hide
p # hide
```

## Related packages

* [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)

    Contains a couple robust statistics, but has no overlapping functionality with this package.

* [RobustStats.jl](https://github.com/mrxiaohe/RobustStats.jl)

    Contains many more robust statistics primarily based on the WRS R package. Appears to be unmaintained and not updated to Julia v1. The `bivar` function is the same as this package's [`midvar`](@ref), although `bivar` does not have definitions for the statistics across axes of an array.

* [astropy.stats](https://github.com/astropy/astropy)

    Python implementations of all the statistics presented here. There are some slight differences in the function signatures, including the default cutoff value `c` (for some statistics).

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/mileslucas/BiweightStats.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/mileslucas/BiweightStats.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/mileslucas/BiweightStats.jl/issues).
