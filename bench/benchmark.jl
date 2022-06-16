using BenchmarkTools
using BiweightStats
using CSV
using DataFrames
using ProgressMeter
using PythonCall
using Random

rng = Random.seed!(110011)

aps = pyimport("astropy.stats")

# define map between functions
fn_map = Dict(
    location => aps.biweight_location,
    scale => aps.biweight_scale,
    midvar => aps.biweight_midvariance,
    midcov => aps.biweight_midcovariance,
    midcor => aps.biweight_midcorrelation
)

fn_names = Dict(
    location => "location",
    scale => "scale",
    midvar => "midvariance",
    midcov => "midcovariance",
    midcor => "midcorrelation"
)
vec_fn = (location, scale, midvar)
cov_fn = (midcov, midcor)


# vector sizes
N = 10 .^ (1:6)

rows = []
@info "Starting benchmark"
@showprogress "Benchmarking" for n in N
    X = randn(rng, n)
    Y = randn(rng, n)
    Z = randn(rng, n)
    D = [X Y Z]

    @showprogress for fn in vec_fn

        t_jl = @belapsed $fn($X)
        py_fn = fn_map[fn]
        t_py = @belapsed $py_fn($X)
        name=fn_names[fn]
        push!(rows, (;n, name, t_jl, t_py))
        @info "$name - $n" t_jl t_py
    end
    @showprogress for fn in cov_fn
        t_jl = @belapsed $fn($X, $Y)
        py_fn = fn_map[fn]
        t_py = @belapsed $py_fn($X, $Y)
        name = fn_names[fn]
        push!(rows, (;n, name, t_jl, t_py))
        @info "$name - $n" t_jl t_py
    end
end

outname =joinpath(@__DIR__, "benchmark_results.csv")
DataFrame(rows) |> CSV.write(outname)
@info "Benchmark complete! Results saved to $outname"
