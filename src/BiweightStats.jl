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
\\forall i \\quad\\mathrm{where}\\quad u_i^2 \\le 1
```

The cutoff factor, ``c``, can be directly related to a Gaussian standard-deviation by multiplying by 1.4826[^2]. So a typical value of ``c=9`` means outliers further than ``13.3\\sigma`` are clipped (for residuals which are truly Gaussian-distributed). In addition, in `BiweightStats`, we also skip `NaN`s and `Inf`s (but not `missing` or `nothing`).

# References

[^1]: [NIST: biweight](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/biweight.pdf)
[^2]: [Median absolute deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation)
# Methods

  - [`location`](@ref)
  - [`scale`](@ref)
  - [`midvar`](@ref)
  - [`midcov`](@ref)
  - [`midcor`](@ref)
"""
module BiweightStats

using Statistics

export BiweightTransform, location, scale, midvar, midcov, midcor

struct BiweightTransform{D,M,C}
    data::D # data filtered to remove NaN
    med::M # median
    cutoff::C # c * MAD
end

"""
    BiweightTransform(X; c=9, M=nothing)

Creates an iterator based on the biweight transform.[^1] This iterator will first filter all input data so that only finite values remain. Then, the iteration will progress using a custom state, which includes a flag to indicate whether the value is within the cutoff, which is `c` times the median-absolute-deviation (MAD). The MAD is based on the deviation from `M`, which will default to the median of `X` if `M` is `nothing`.

!!! note "Advanced usage"
    
    This transform iterator is used for the internal calculations in `BiweightStats.jl`, which is why it has a somewhat complicated iterator implementation.

# Examples

```jldoctest transform
julia> X = randn(rng, 100);


julia> X[10] = 1e4 # add clear outlier
10000.0

julia> X[13] = NaN # add NaN
NaN

julia> X[25] = Inf # add Inf
Inf

julia> bt = BiweightTransform(X);

```

Lets confirm all the entries are finite. The iteration interface is divided into

```julia
(d, u2, flag), state = iterate(bt, [state])
```

where `d` is the data value minus `M`, `u2` is `(d / (c * MAD))^2`, and `flag` is whether the value is within the transformed dataset.

```jldoctest transform
julia> all(d -> isfinite(d[1]), bt)
true
```

and let's see how iteration differs between a normal sample and an outlier sample, which we manually inserted at index `10`-

```jldoctest transform
julia> (d, u2, flag), _ = iterate(bt, 9)
((-0.17093842061187192, 0.0009098761083851183, true), 10)

julia> (d, u2, flag), _ = iterate(bt, 10)
((0.0, 0.0, false), 11)
```

# References

[^1]: [NIST: biweight](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/biweight.pdf)
"""
function BiweightTransform(X; c=9, M=nothing)
    data = filter(isfinite, X)
    if isnothing(M)
        med = median(data)
    else
        med = M
    end
    mad = median!(abs.(data .- med))
    return BiweightTransform(data, med, c * mad)
end

Base.size(bt::BiweightTransform) = size(bt.data)
Base.size(bt::BiweightTransform, d) = size(bt.data, d)
Base.IteratorSize(bt::BiweightTransform) = Base.IteratorSize(bt.data)
Base.IteratorEltype(bt::BiweightTransform) = Base.IteratorEltype(bt.data)

function Base.iterate(bt::BiweightTransform)
    next = iterate(bt.data)
    isnothing(next) && return nothing
    data, state = next
    return _biweight_iterate(bt, data), state
end

function Base.iterate(bt::BiweightTransform, state)
    next = iterate(bt.data, state)
    isnothing(next) && return nothing
    data, state = next
    return _biweight_iterate(bt, data), state
end

function _biweight_iterate(bt, data)
    d = data - bt.med
    u2 = (d / bt.cutoff)^2 # ((x - Mx) / (c * MAD))^2
    # if MAD is zero, u2 will become 0
    iszero(bt.cutoff) && return zero(d), zero(u2), true
    # filter values beyond cutoff
    if u2 > 1
        return zero(d), zero(u2), false
    end
    # normal return
    return d, u2, true
end

"""
    location(X; c=9, M=nothing)
    location(X::AbstractArray; dims=:, kwargs...)

Calculate the biweight location, a robust measure of location.

```math
\\hat{y} = \\frac{\\sum_{u_i^2 \\le 1}{y_i(1 - u_i^2)^2}}{\\sum_{u_i^2 \\le 1}{(1 - u_i^2)^2}}
```

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000) .+ 50;


julia> location(X)
49.98008021468018
```

## Iterative refinement

You can iteratively refine the location estimate by manually passing the median, like so-

```jldoctest
X = 10 .* randn(rng, 1000) .+ 50
let ystar, ystar_old
    ystar = ystar_old = location(X)
    tol = 1e-6
    maxiter = 10
    for _ = 1:maxiter
        ystar = location(X; M=ystar_old)
        isapprox(ystar_old, ystar; atol=tol) && break
        ystar_old = ystar
    end
    ystar
end

# output

49.991666155308145
```

# References

 1. [NIST: biweight location](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwloc.htm)
"""
location(X; maxiter=10, tol=1e-6, kwargs...) = biweight_location(X; maxiter, tol, kwargs...)
location(X::AbstractArray; dims=:, kwargs...) = location(X, dims; kwargs...)
location(X::AbstractArray, ::Colon; kwargs...) = biweight_location(X; kwargs...)
function location(X::AbstractArray, dims::Int; kwargs...)
    return mapslices(sl -> biweight_location(sl; kwargs...), X; dims)
end

function biweight_location(X; kwargs...)
    T = float(eltype(X))
    itr = BiweightTransform(X; kwargs...)
    num = zero(T)
    den = zero(T)
    for (d, u2, flag) in itr
        flag || continue
        w = (1 - u2)^2
        num += w * (d + itr.med)
        den += w
    end
    return num / den
end

"""
    scale(X; c=9, M=nothing)
    scale(X::AbstractArray; dims=:, kwargs...)

Compute the biweight scale of the variable. This is the same as the square-root of the midvariance.

```math
\\hat{\\sigma} = \\frac{\\sqrt{n\\sum^n_{u_i^2 \\le 1}{(y_i - \\bar{y})^2(1 - u_i^2)^4}}}{\\sum^n_{u_i^2 \\le 1}{(1 - u_i^2)(1 - 5u_i^2)}}
```

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000) .+ 50;


julia> scale(X)
10.045813567765071
```

# References

 1. [NIST: biweight scale](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwscale.htm)

# See Also

[`midcor`](@ref), [`midvar`](@ref), [`midcov`](@ref)
"""
scale(X; kwargs...) = sqrt(biweight_midvar(X; kwargs...))
scale(X::AbstractArray; dims=:, kwargs...) = scale(X, dims; kwargs...)
scale(X::AbstractArray, ::Colon; kwargs...) = sqrt(biweight_midvar(X; kwargs...))
function scale(X::AbstractArray, dims::Int; kwargs...)
    return mapslices(sl -> sqrt(biweight_midvar(sl; kwargs...)), X; dims)
end

"""
    midvar(X; c=9, M=nothing)
    midvar(X::AbstractArray; dims=:, kwargs...)

Compute the biweight midvariance of the variable.

```math
\\hat{\\sigma^2} = \\frac{n\\sum^n_{u_i^2 \\le 1}{(y_i - \\bar{y})^2(1 - u_i^2)^4}}{\\left[\\sum_{u_i^2 \\le 1}{(1 - u_i^2)(1 - 5u_i^2)}\\right]^2}
```

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000) .+ 50;


julia> midvar(X)
100.9183702382928
```

# References

 1. [NIST: biweight midvariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm)

# See Also

[`scale`](@ref), [`midcor`](@ref), [`midcov`](@ref)
"""
midvar(X; kwargs...) = biweight_midvar(X; kwargs...)
midvar(X::AbstractArray; dims=:, kwargs...) = midvar(X, dims; kwargs...)
midvar(X::AbstractArray, ::Colon; kwargs...) = biweight_midvar(X; kwargs...)
function midvar(X::AbstractArray, dims::Int; kwargs...)
    return mapslices(sl -> biweight_midvar(sl; kwargs...), X; dims)
end

function biweight_midvar(X; kwargs...)
    itr = BiweightTransform(X; kwargs...)
    # init
    T = eltype(X)
    num = zero(T)
    den = zero(T)
    n = 0
    for (d, u2, flag) in itr
        flag || continue
        num += d^2 * (1 - u2)^4
        den += (1 - u2) * (1 - 5 * u2)
        n += 1
    end

    return n * num / den^2
end

"""
    midcov(X, [Y]; c=9)

Computes biweight midcovariance between the two vectors. If only one vector is provided the biweight midvariance will be calculated.

```math
\\hat{\\sigma}_{xy} = \\frac{\\sum_{u_i^2 \\le 1,v_i^2 \\le 1}{(x_i - \\bar{x})(1 - u_i^2)^2(y_i - \\bar{y})(1 - v_i^2)^2}}{\\sum_{u_i^2 \\le 1}{(1 - u_i^2)(1 - 5u_i^2)}\\sum_{v_i^2 \\le 1}{(1 - v_i^2)(1 - 5v_i^2)}}
```

!!! warning
    
    `NaN` and `Inf` cannot be removed in the covariance calculation, so if they are present the returned value will be `NaN`. To prevent this, consider imputing values for the non-finite data.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000, 2) .+ 50;


julia> midcov(X[:, 1], X[:, 2])
-1.058463590812247

julia> midcov(X[:, 1]) ≈ midvar(X[:, 1])
true

julia> X[3, 2] = NaN;


julia> midcov(X[:, 1], X[:, 2])
NaN
```

# References

 1. [NIST: biweight midcovariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm)

# See Also

[`scale`](@ref), [`midvar`](@ref), [`midcor`](@ref)
"""
function midcov(X, Y; kwargs...)
    length(X) == length(Y) || throw(DimensionMismatch("collections must have equal length"))
    T = float(promote_type(eltype(X), eltype(Y)))
    (all(isfinite, X) && all(isfinite, Y)) || return T(NaN)
    itrx = BiweightTransform(X; kwargs...)
    itry = BiweightTransform(Y; kwargs...)
    # init
    num = zero(T)
    den1 = zero(T)
    den2 = zero(T)
    n = 0
    for ((dx, u2, flagx), (dy, v2, flagy)) in zip(itrx, itry)
        (flagx && flagy) || continue
        num += dx * (1 - u2)^2 * dy * (1 - v2)^2
        den1 += (1 - u2) * (1 - 5 * u2)
        den2 += (1 - v2) * (1 - 5 * v2)
        n += 1
    end

    return n * num / (den1 * den2)
end

midcov(X; kwargs...) = midvar(X; kwargs...)

"""
    midcov(X::AbstractMatrix; dims=1, c=9)

Computes the variance-covariance matrix using the biweight midcovariance. By default, each column is a separate variable, so an `(M, N)` matrix with `dims=1` will create an `(N, N)` covariance matrix. If `dims=2`, though, each row will become a variable, leading to an `(M, M)` covariance matrix.

!!! warning
    
    `NaN` and `Inf` cannot be removed in the covariance calculation, so if they are present the returned value will be `NaN`. To prevent this, consider imputing values for the non-finite data.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000, 3) .+ 50;


julia> C = midcov(X)
3×3 Matrix{Float64}:
 100.918    -1.05846    -2.88515
  -1.05846  94.702      -0.490742
  -2.88515  -0.490742  100.699

julia> size(midcov(X; dims=2))
(1000, 1000)
```

# References

 1. [NIST: biweight midcovariance](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm)

# See Also

[`scale`](@ref), [`midvar`](@ref), [`midcor`](@ref)
"""
function midcov(X::AbstractMatrix{T}; dims=1, kwargs...) where {T}
    vardim = dims == 1 ? 2 : 1
    out = zeros(float(T), size(X, vardim), size(X, vardim))

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
    midcor(X, Y; c=9)

Compute the correlation between two variables using the midvariance and midcovariances.

```math
\\frac{s_{xy}}{\\sqrt{s_{xx} \\cdot s_{yy}}}
```

where ``s_{xx},s_{yy}`` are the midvariances of each vector, and ``s_{xy}`` is the midcovariance of the two vectors.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000, 2) .+ 50;


julia> midcor(X[:, 1], X[:, 2])
-0.010827077678217934
```

# References

 1. [Wikipedia](https://en.wikipedia.org/wiki/Biweight_midcorrelation)
 2. [NIST: Biweight midcorrelation](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidcr.htm)

# See Also

[`midvar`](@ref), [`midcov`](@ref), [`scale`](@ref)
"""
function midcor(X, Y; kwargs...)
    sxy = midcov(X, Y; kwargs...)
    sxx = midcov(X; kwargs...)
    syy = midcov(Y; kwargs...)
    return sxy / sqrt(sxx * syy)
end

"""
    midcor(X::AbstractMatrix; dims=1, c=9)

Computes the correlation matrix using the biweight midcorrealtion. By default, each column of the matrix is a separate variable, so an `(M, N)` matrix with `dims=1` will create an `(N, N)` correlation matrix. If `dims=2`, though, each row will become a variable, leading to an `(M, M)` correlation matrix. The diagonal will always be one.

# Examples

```jldoctest
julia> X = 10 .* randn(rng, 1000, 3) .+ 50;


julia> C = midcor(X)
3×3 Matrix{Float64}:
  1.0        -0.0108271  -0.0286201
 -0.0108271   1.0        -0.0050253
 -0.0286201  -0.0050253   1.0

julia> size(midcor(X; dims=2))
(1000, 1000)
```

# References

 1. [Wikipedia](https://en.wikipedia.org/wiki/Biweight_midcorrelation)
 2. [NIST: Biweight midcorrelation](https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidcr.htm)

# See Also

[`midvar`](@ref), [`midcov`](@ref), [`scale`](@ref)
"""
function midcor(X::AbstractMatrix{V}; dims=1, kwargs...) where {V}
    vardim = dims == 1 ? 2 : 1
    T = float(V)
    out = zeros(T, size(X, vardim), size(X, vardim))

    for i in axes(out, 1), j in axes(out, 2)
        if i > j
            out[i, j] = out[j, i]
            continue
        end
        if i == j
            cor = one(T)
        else
            x = selectdim(X, vardim, i)
            y = selectdim(X, vardim, j)
            cor = midcor(x, y; kwargs...)
        end
        out[i, j] = cor
    end

    return out
end

end
