# BiweightStats.jl

[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BiweightStats.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Build Status](https://github.com/mileslucas/BiweightStats.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mileslucas/BiweightStats.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mileslucas/BiweightStats.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mileslucas/BiweightStats.jl)
[![License](https://img.shields.io/github/license/mileslucas/BiweightStats.jl?color=yellow)](LICENSE)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mileslucas.github.io/BiweightStats.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mileslucas.github.io/BiweightStats.jl/dev)

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

## Methods

The following statistics are defined as methods. For more information, see the [online documentation](https://mileslucas.github.io/BiweightStats.jl/dev)

* `location`
* `scale`
* `midvar`
* `midcov`
* `midcor`

## Related packages

* [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)

    Contains a couple robust statistics, but has no overlapping functionality with this package.

* [RobustStats.jl](https://github.com/mrxiaohe/RobustStats.jl)

    Contains many more robust statistics primarily based on the WRS R package. Appears to be unmaintained and not updated to Julia v1. The `bivar` function is the same as this package's [`midvar`](@ref), although `bivar` does not have definitions for the statistics across axes of an array.

* [astropy.stats](https://github.com/astropy/astropy)

    Python implementations of all the statistics presented here. Some slight differences in the function signatures and the implementations are independent.

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/mileslucas/BiweightStats.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/mileslucas/BiweightStats.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/mileslucas/BiweightStats.jl/issues).
