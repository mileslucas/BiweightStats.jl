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

using Statistics

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
function midvariance(values; c=6)
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


end
