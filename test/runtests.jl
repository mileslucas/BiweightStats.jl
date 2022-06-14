using BiweightStats
using StableRNGs
using Test

rng = StableRNG(1123)

@testset "BiweightStats.jl" begin
    @testset "BiweightTransform" begin
        X = randn(rng, 100)
        X[12] = 1e5
        X[27] = NaN
        bt = BiweightStats.BiweightTransform(X)
        @test all(isfinite, bt.data)
        (d, u2, flag), state = iterate(bt, 12)
        @test state == 13
        @test !flag
        @test iszero(u2) && iszero(d)
    end

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
        @test val ≈ 2 atol=1e-5
    end

    @testset "scale" begin
        X = randn(rng, 10000)
        val = scale(X)
        # depend on midvar tests
        @test val == sqrt(midvar(X))
        # axes
        # X = randn(rng, 10000, 5)
        # for dims in (1, 2)
        #     vals = scale(X; dims=dims)
        #     @test vals ≈ sqrt.(midvar(X, dims=dims))
        # en
    end

    @testset "midvar" begin
        X = randn(rng, 10000)
        val = midvar(X)
        @test val ≈ 1 atol=0.03

        # trivial variance
        @test midvar(zeros(100)) ≈ midvar(ones(100)) ≈ 0

        # axes
        # X = randn(rng, 10000, 5)
        # vals = midvar(X; dims=1)
        # @test vals ≈ [1 1 1 1] rtol=1e-2
        # vals = midvar(X; dims=2)
        # @test size(vals) == (10000, 1)

        # outlier
        val = midvar([1, 2, 3, 500, 2])
        @test val ≈ 0.51266 atol=1e-5

        # NaN
        val = midvar([1, 2, 3, NaN, 2])
        @test val ≈ 0.55472 atol=1e-5
    end
end
