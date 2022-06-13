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

using DocStringExtensions
using Statistics

"""
    $(SIGNATURES)


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
        isnan(d) && continue
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
function midvar(values; c=9)
    Δ = values .- median(values)
    mad = median(Iterators.map(abs, Δ))
    # init
    T = eltype(values)
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
$(SIGNATURES)
"""

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

function midcov(X::AbstractVector{V}, Y::AbstractVector{S}; c=9) where {V,S}
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
$(SIGNATURES)
"""
function midcor(X, Y=X; kwargs...)
    cov = midcov(X, Y; kwargs...)
    return cov[begin, end] / sqrt(cov[begin, begin] * cov[end, end])
end

end
