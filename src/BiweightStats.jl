"""
    BiweightStats

A module for robust statistics based on the biweight transform.[^1]

# Biweight Statistics

The basis of the biweight transform is robust analysis, that is, statistics which are resilient to outliers while still efficiently representing a variety of underlying distributions. The biweight transform is based off the *median* and the *median absolute deviation (MAD)*. The median is a robust estimator of location, and the MAD is a robust estimator of scale

```math
\\mathrm{MAD}(X) = \\mathrm{median}\\left|X_i - \\bar{X}\\right|
```
where ``\\bar{X}`` is the median.

The biweight transform improves upon these estimates by filtering out data beyond a critical cutoff. The analogy is doing a sigma-filter, but using these robust statistics instead of the standard deviation and mean.

```math
u_i = \\frac{X_i - \\bar{X}}{c \\cdot \\mathrm{MAD}}
```

```math
\\forall i \\quad\\mathrm{where}\\quad u_i^2 < 1
```

The cutoff factor, ``c``, can be directly related to a Gaussian standard-deviation by multiplying by 1.48. So a typical value of ``c=6`` means outliers further than ``\\sim 9\\sigma`` are clipped. In addition, in `BiweightStats`, we also skip `NaN`s.

# References

[^1]: [NIST: biweight](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/biweight.pdf)

# Methods

- [`location`](@ref)
- [`scale`](@ref)
- [`midvar`](@ref)
- [`midcov`](@ref)
- [`midcor`](@ref)
"""
module BiweightStats

using Statistics

export location, scale, midvar, midcov, midcor

"""
    location(X; c=6, maxiter=10, tol=1e-6)

Iteratively calculate the biweight location, a robust measure of location.

# Stopping Criteria

The location will be refined until `maxiter` is reached or until the absolute change between estimates is below `tol`.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 10) .+ 50;

julia> location(X)
49.483640367304304
```

# References

1. [NIST: biweight location](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwloc.htm)
"""
function location(X; c=6, maxiter=10, tol=1e-6)
    T = eltype(X)
    ystar = zero(T)
    ystar_old = ystar
    num = zero(T)
    den = zero(T)
    for _ in 1:maxiter
        Δ = X .- ystar
        mad = median(Iterators.map(abs, Iterators.filter(!isnan, Δ)))
        num = zero(T)
        den = zero(T)
        for (y, d) in zip(X, Δ)
            isnan(d) && continue
            u2 = (d / (c * mad))^2
            u2 > 1 && continue
            w = (1 - u2)^2
            num += w * y
            den += w
        end
        ystar_old = ystar
        ystar = num / den
        isapprox(ystar, ystar_old, atol=tol) && break
    end
    return ystar
end

"""
    scale(X; c=6)

Compute the biweight midvariance of the variable.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 10) .+ 50;

julia> scale(X)
9.269099543377779
```

# References

1. [NIST: biweight scale](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwscale.htm)

# See Also

[`midcor`](@ref), [`midvar`](@ref), [`midcov`](@ref)
"""
function scale(X; c=6)
    Δ = X .- median(X)
    mad = median(Iterators.map(abs, Iterators.filter(!isnan, Δ)))
    # init
    T = eltype(X)
    num = zero(T)
    den = zero(T)
    n = 0
    for d in Δ
        isnan(d) && continue
        u2 = (d / (c * mad))^2
        u2 > 1 && continue
        num += d^2 * (1 - u2)^4
        den += (1 - u2) * (1 - 5 * u2)
        n += 1
    end

    return sqrt(n * num / (den * (den - 1)))
end

"""
    midvar(X; c=6)

Compute the biweight midvariance of the variable.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 10) .+ 50;

julia> midvar(X)
70.53324162560791
```

# References

1.  [NIST: biweight midvariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm)

# See Also

[`scale`](@ref), [`midcor`](@ref), [`midcov`](@ref)
"""
function midvar(X; c=6)
    Δ = X .- median(X)
    mad = median(Iterators.map(abs, Iterators.filter(!isnan, Δ)))
    # init
    T = eltype(X)
    num = zero(T)
    den = zero(T)
    n = 0
    for d in Δ
        isnan(d) && continue
        u2 = (d / (c * mad))^2
        u2 > 1 && continue
        num += d^2 * (1 - u2)^4
        den += (1 - u2) * (1 - 5 * u2)
        n += 1
    end

    return n * num / den^2
end

"""
    midcov(X::AbstractVector, [Y::AbstractVector]; c=6)

Computes biweight midcovariance between the two vectors. If only one vector is provided the biweight midvariance will be calculated.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 10, 2) .+ 50;

julia> midcov(X[:, 1], X[:, 2])
-25.644749715907157

julia> midcov(X[:, 1]) ≈ midvar(X[:, 1])
true
```

# References

1. [NIST: biweight midcovariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm)

# See Also

[`scale`](@ref), [`midvar`](@ref), [`midcor`](@ref)
"""
function midcov(X::AbstractVector{V}, Y::AbstractVector{S}; c=6) where {V,S}
    length(X) == length(Y) || throw(DimensionMismatch("vectors must have same length"))
    Δx = X .- median(X)
    Δy = Y .- median(Y)
    madx = median(Iterators.map(abs, Iterators.filter(!isnan, Δx)))
    mady = median(Iterators.map(abs, Iterators.filter(!isnan, Δy)))
    # init
    T = promote_type(V, S)
    num = zero(T)
    den1 = zero(T)
    den2 = zero(T)
    n = 0
    for (dx, dy) in zip(Δx, Δy)
        (isnan(dx) || isnan(dy)) && continue
        u2 = (dx / (c * madx))^2
        v2 = (dy / (c * mady))^2
        if u2 ≤ 1
            den1 += (1 - u2) * (1 - 5 * u2)
        end
        if v2 ≤ 1
            den2 += (1 - v2) * (1 - 5 * v2)
        end
        if u2 ≤ 1 && v2 ≤ 1
            num += dx * (1 - u2)^2 * dy * (1 - v2)^2
            n += 1
        end
    end

    return n * num / (den1 * den2)
end

midcov(X; kwargs...) = midvar(X; kwargs...)

"""
    midcov(X::AbstractMatrix; dims=1, c=6)

Computes the variance-covariance matrix using the biweight midcovariance. By default, each column is a separate variable, so an `(M, N)` matrix with `dims=1` will create an `(N, N)` covariance matrix. If `dims=2`, though, each row will become a variable, leading to an `(M, M)` covariance matrix.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 5, 3) .+ 50;

julia> C = midcov(X)
3×3 Matrix{Float64}:
 214.707   -38.4913   38.2848
 -38.4913   35.6904  -18.5234
  38.2848  -18.5234   77.4838

julia> midcov(X; dims=2)
5×5 Matrix{Float64}:
  7.10619   13.7261     26.8059    2.28917   -0.640045
 13.7261    26.6125     52.0541    1.31522   -0.367732
 26.8059    52.0541    327.443   -83.8513    70.6212
  2.28917    1.31522   -83.8513  128.887    -53.5744
 -0.640045  -0.367732   70.6212  -53.5744    29.7501
```

# References

1. [NIST: biweight midcovariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm)

# See Also

[`scale`](@ref), [`midvar`](@ref), [`midcor`](@ref)
"""
function midcov(X::AbstractMatrix{T}; dims=1, kwargs...) where T
    vardim = dims == 1 ? 2 : 1
    out = zeros(T, size(X, vardim), size(X, vardim))

    for i in axes(out, 1), j in axes(out, 2)
        if i > j
            out[i, j] = out[j, i]
            continue
        end
        if i == j
            x = selectdim(X, vardim, i)
            cov = midcov(x; kwargs...)
        else
            x = selectdim(X, vardim, i)
            y = selectdim(X, vardim, j)
            cov = midcov(x, y; kwargs...)
        end
        out[i, j] = cov
    end

    return out
end

"""
    midcor(X::AbstractVector, Y::AbstractVector; c=6)

Compute the correlation between two variables using the midvariance and midcovariances.

```math
\\frac{s_{xy}}{\\sqrt{s_{xx} \\cdot s_{yy}}}
```
where ``s_{xx},s_{yy}`` are the midvariances of each vector, and ``s_{xy}`` is the midcovariance of the two vectors.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 10, 2) .+ 50;

julia> midcor(X[:, 1], X[:, 2])
-0.22014156692229941
```

# References

1. [Wikipedia](https://en.wikipedia.org/wiki/Biweight_midcorrelation)
2. [NIST: Biweight midcorrelation](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidcr.htm)

# See Also

[`midvar`](@ref), [`midcov`](@ref), [`scale`](@ref)
"""
function midcor(X::AbstractVector, Y::AbstractVector; kwargs...)
    sxx = midcov(X; kwargs...)
    syy = midcov(Y; kwargs...)
    sxy = midcov(X, Y; kwargs...)
    return sxy / sqrt(sxx * syy)
end

end
