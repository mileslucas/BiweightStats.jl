# Benchmarks

Here you can find benchmarks of the BiweightStats.jl estimators versus the [astropy](https://github.com/astropy/astropy) implementations in `astropy.stats`.

## Setup

The benchmark is written in Julia ([./benchmark.jl]()) and uses [PythonCall.jl](https://github.com/cjdoris/PythonCall.jl) for interfacing with `astropy.stats. All you need to do is activate and instantiate the Julia environment, and the package managers should handle the rest

```julia
julia>] activate bench

(bench) pkg> instantiate
```

**Note**: You only need to instantiate the environment if it's your first time running the code.

## Usage

From the Julia REPL (make sure you've activated the `bench` environment)

```julia
julia> include("bench/benchmark.jl")
```

or from the command-line directly

```
julia --project=bench bench/benchmark.jl
```

## Results

The results are saved in a CSV file ([./benchmark_results.csv]()) and are plotted in the [online documentation](https://mileslucas.github.io/BiweightStats.jl/dev/index.html#Benchmarks)
