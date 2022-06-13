"""
    Biweight

A module for robust statistics. The following methods are implemented

# Methods

- [`Biweight.location`](@ref)
- [`Biweight.scale`](@ref)
- [`Biweight.midvar`](@ref)
- [`Biweight.midcov`](@ref)
- [`Biweight.midcor`](@ref)
"""
module Biweight

using Statistics

"""
    Biweight.location(X; c=6, maxiter=10, tol=1e-6)

Iteratively calculate the biweight location, a robust measure of location.

# Stopping Criteria

The location will be refined until `maxiter` is reached or until the absolute change between estimates is below `tol`.

"""
function location(X; c=6, maxiter=10, tol=1e-6)
    T = eltype(X)
    ystar = zero(T)
    ystar_old = ystar
    num = zero(T)
    den = zero(T)
    for _ in 1:maxiter
        Δ = X .- ystar
        mad = median(Iterators.map(abs, Δ))
        num = zero(T)
        den = zero(T)
        for (y, d) in zip(X, Δ)
            isnan(d) && continue
            u2 = (d / (c * mad))^2
            u2 ≥ 1 && continue
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
    Biweight.scale(X; c=9)

Compute the biweight midvariance of the variable.

# References

1. https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm

# See Also

[`Biweight.midcor`](@ref), [`Biweight.midvar`](@ref), [`Biweight.midcov`](@ref)
"""
function scale(X; c=9)
    Δ = X .- median(X)
    mad = median(Iterators.map(abs, Δ))
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
    Biweight.midvar(X; c=9)

Compute the biweight midvariance of the variable.

# References

1. https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm

# See Also

[`Biweight.scale`](@ref), [`Biweight.midcor`](@ref), [`Biweight.midcov`](@ref)
"""
function midvar(X; c=9)
    Δ = X .- median(X)
    mad = median(Iterators.map(abs, Δ))
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

function midcov(X::AbstractVector{T}; c=9) where T
    Δ = X .- median(X)
    mad = median(Iterators.map(abs, Δ))
    # init
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
    Biweight.midcov(X::AbstractVecotr, [Y::AbstractVector]; c=9)

Computes biweight midcovariance between the two vectors. If only one vector is provided the biweight midvariance will be calculated.

# References

1. https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm

# See Also

[`Biweight.scale`](@ref), [`Biweight.midvar`](@ref), [`Biweight.midcor`](@ref)
"""
function midcov(X::AbstractVector{V}, Y::AbstractVector{S}; c=9) where {V,S}
    length(X) == length(Y) || throw(DimensionMismatch("vectors must have same length"))
    Δx = X .- median(X)
    Δy = Y .- median(Y)
    madx = median(Iterators.map(abs, Δx))
    mady = median(Iterators.map(abs, Δy))
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
        if u2 < 1
            den1 += (1 - u2) * (1 - 5 * u2)
        end
        if v2 < 1
            den2 += (1 - v2) * (1 - 5 * v2)
        end
        if u2 < 1 && v2 < 1
            num += dx * (1 - u2)^2 * dy * (1 - v2)^2
            n += 1
        end
    end

    return n * num / (den1 * den2)
end


"""
    Biweight.midcov(X::AbstractMatrix; dims=1, c=9)

Computes the variance-covariance matrix using the biweight midcovariance. By default, each column is a separate variable, so an `(M, N)` matrix with `dims=1` will create an `(N, N)` covariance matrix. If `dims=2`, though, each row will become a variable, leading to an `(M, M)` covariance matrix.

# References

1. https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm

# See Also

[`Biweight.scale`](@ref), [`Biweight.midvar`](@ref), [`Biweight.midcor`](@ref)
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

raw"""
    Biweight.midcor(X::AbstractVector, Y::AbstractVector; c=9)

Compute the correlation between two variables using the midvariance and midcovariances.

```math
\frac{\zeta_{xy}}{\zeta_{xx}\zeta_{yy}}
```
where $\zeta_{xx},\zeta_{yy}$ are the midvariances of each vector, and $\zeta_{xy}$ is the midcovariance of the two vectors.

# References

1. https://en.wikipedia.org/wiki/Biweight_midcorrelation
2. https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidcr.htm

# See Also

[`Biweight.midvar`](@ref), [`Biweight.midcov`](@ref), [`Biweight.scale`](@ref)
"""
function midcor(X::AbstractVector, Y::AbstractVector; kwargs...)
    sxx = midcov(X; kwargs...)
    syy = midcov(X; kwargs...)
    sxy = midcov(X, Y; kwargs...)
    return sxy / sqrt(sxx * syy)
end

end
