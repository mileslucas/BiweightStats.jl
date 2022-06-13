```@meta
CurrentModule = Biweight
```

# Biweight.jl

[![Code](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/mileslucas/Biweight.jl)
[![Build Status](https://github.com/mileslucas/Biweight.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mileslucas/Biweight.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/Biweight.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/mileslucas/Biweight.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mileslucas/Biweight.jl)
[![License](https://img.shields.io/github/license/mileslucas/Biweight.jl?color=yellow)](https://github.com/mileslucas/Biweight.jl/blob/main/LICENSE)

Robust statistics based on the biweight transform.

## Installation

Biweight.jl can be installed using the Julia package manager. From the Julia REPL, enter Pkg mode (by pressing `]`)

```julia
julia>]

pkg> add https://github.com/mileslucas/Biweight.jl
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with

```julia
julia> using Biweight
```

To exit Pkg mode, just backspace. Once the package is installed it can be imported with
For more information, see the [Pkg documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Related packages

1. [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)

    Contains a couple robust statistics, but has no overlapping functionality with this package.

2. [RobustStats.jl](https://github.com/mrxiaohe/RobustStats.jl)

    Contains many more robust statistics primarily based on the WRS R package. Appears to be unmaintained and not updated to Julia v1. The `bivar` function is the same as this package's [`Biweight.midvar`](@ref), although `bivar` does not have definitions for the statistics across axes of an array.

## API/Reference

```@index
```

```@docs
Biweight.location
Biweight.scale
Biweight.midvar
Biweight.midcov
Biweight.midcor
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/mileslucas/Biweight.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/mileslucas/Biweight.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/mileslucas/Biweight.jl/issues).
