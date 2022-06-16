using BenchmarkTools
using BiweightStats
using CSV
using ProgressLogging
using PythonCall
using Random

rng = Random.seed!(110011)

aps = pyimport("astropy.stats")

# vector sizes
N = 10 .^ (1:5)

rows = []
@info "Starting benchmark"
@progress "Benchmarking" for n in N
    X = randn(rng, n)
    Y = randn(rng, n)

    t_jl = @belapsed location($X; c=9)
    t_py = @belapsed aps.biweight_location($X; c=9)
    push!(rows, (;n, name="location", t_jl, t_py))
    @info "location - $n" t_jl t_py

    t_jl = @belapsed scale($X; c=9)
    t_py = @belapsed aps.biweight_scale($X; c=9, modify_sample_size=true)
    push!(rows, (;n, name="scale", t_jl, t_py))
    @info "scale - $n" t_jl t_py

    t_jl = @belapsed midvar($X; c=9)
    t_py = @belapsed aps.biweight_midvariance($X; c=9, modify_sample_size=true)
    push!(rows, (;n, name="midvariance", t_jl, t_py))
    @info "midvariance - $n" t_jl t_py

    Z = [X Y]
    t_jl = @belapsed midcov($Z; c=9)
    t_py = @belapsed aps.biweight_midcovariance($(Z'); c=9, modify_sample_size=true)
    push!(rows, (;n, name="midcovariance", t_jl, t_py))
    @info "midcovariance - $n" t_jl t_py

    t_jl = @belapsed midcor($X, $Y; c=9)
    t_py = @belapsed aps.biweight_midcorrelation($X, $Y; c=9, modify_sample_size=true)
    push!(rows, (;n, name="midcorrelation", t_jl, t_py))
    @info "midcorrelation - $n" t_jl t_py
end

outname =joinpath(@__DIR__, "benchmark_results.csv")
CSV.write(outname, rows)
@info "Benchmark complete! Results saved to $outname"
