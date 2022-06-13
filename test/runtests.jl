using BiweightStats
using StableRNGs
using Test

rng = StableRNG(1123)

@testset "BiweightStats.jl" begin
    @testset "location" begin
        X = randn(rng, 10000)
        val = location(X)
        @test val ≈ 0 atol=1e-2

        X = ones(10)
        val = location(X)
        @test val ≈ 1

        # X = repeat([1 2 3], outer=5)
        # vals = location(X; dims=1)
        # @test vals ≈ [1 2 3]
        # vals = location(X; dims=2)
        # @test vals ≈ fill(2.0, 5)

        # outlier
        val = location([1, 2, 3, 500, 2])
        @test val ≈ 2

        # NaN
        val = location([1, 2, 3, NaN, 2])
        @test val ≈ 2
    end

    @testset "scale" begin
        X = randn(rng, 10000)
        val = scale(X)
        @test val ≈ 1 rtol=1e-2
        # @test val ≈ sqrt(midvar(X))

        # X = randn(rng, 10000, 5)
        # vals = scale(X; dims=1)
        # @test vals ≈ [1 1 1 1] rtol=1e-2
        # vals = scale(X; dims=2)
        # @test size(vals) == (10000, 1)

        # outlier
        val = scale([1, 2, 3, 500, 2])
        @test val ≈ 0.85285 atol=1e-5

        # NaN
        val = scale([1, 2, 3, NaN, 2])
        @test val ≈ 0.85285 atol=1e-5
    end

    @testset "midvar" begin
        X = randn(rng, 10000)
        val = midvar(X)
        @test val ≈ 1 rtol=1e-2

        # X = randn(rng, 10000, 5)
        # vals = midvar(X; dims=1)
        # @test vals ≈ [1 1 1 1] rtol=1e-2
        # vals = midvar(X; dims=2)
        # @test size(vals) == (10000, 1)

        # outlier
        val = midvar([1, 2, 3, 500, 2])
        @test val ≈ 0.5294 atol=1e-5

        # NaN
        val = midvar([1, 2, 3, NaN, 2])
        @test val ≈ 0.5294 atol=1e-5
    end
end
