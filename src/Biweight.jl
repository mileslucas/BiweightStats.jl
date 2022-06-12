"""
    Biweight

A module for robust statistics. The following methods are implemented

# Methods

- [`Biweight.location`](@ref)
- [`Biweight.scale`](@ref)
- [`Biweight.midvariance`](@ref)
- [`Biweight.midcovariance`](@ref)
- [`Biweight.midcorrelation`](@ref)
"""
module Biweight

using DocStringExtensions
using Statistics

"""
$(SIGNATURES)
"""
function location(values; c=6, iters=10, tol=1e-6)
    T = eltype(values)
    ystar = zero(T)
    ystar_old = ystar
    num = zero(T)
    den = zero(T)
    for _ in 1:iters
        Δ = values .- ystar
        mad = median(Iterators.map(abs, Δ))
        num = zero(T)
        den = zero(T)
        for (y, d) in zip(values, Δ)
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
$(SIGNATURES)
"""
function scale(values; c=9)
    Δ = values .- median(values)
    mad = median(Iterators.map(abs, Δ))
    # init
    T = eltype(values)
    num = zero(T)
    den = zero(T)
    n = 0
    for d in Δ
        u2 = (d / (c * mad))^2
        u2 > 1 && continue
        num += d^2 * (1 - u2)^4
        den += (1 - u2) * (1 - 5 * u2)
        n += 1
    end

    return sqrt(n * num / (den * (den - 1)))
end

# biweight midvariance
"""
$(SIGNATURES)
"""
function midvariance(values; c=9)
    Δ = values .- median(values)
    mad = median(Iterators.map(abs, Δ))
    # init
    T = eltype(values)
    num = zero(T)
    den = zero(T)
    n = 0
    for d in Δ
        u2 = (d / (c * mad))^2
        u2 > 1 && continue
        num += d^2 * (1 - u2)^4
        den += (1 - u2) * (1 - 5 * u2)
        n += 1
    end

    return n * num / den^2
end

"""
$(SIGNATURES)
"""
function midcovariance(X, Y=X; c=9)
end

"""
$(SIGNATURES)
"""
function midcorrelation(X, Y=X; kwargs...)
    cov = midcovariance(X, Y; kwargs...)
    return cov[begin, end] / sqrt(cov[begin, begin] * cov[end, end])
end

end
